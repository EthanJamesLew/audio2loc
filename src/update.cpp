#include "update.h"

#include <zlib.h>

#include <condition_variable>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <filesystem>
#include <thread>

#include "io.h"
#include "util.h"

namespace fs = std::filesystem;

namespace pl {

std::string utc_date(double t) {
  return fmt_time(std::floor(t / 86400) * 86400).substr(0, 10);
}

std::vector<std::string> dates_between(double t0, double t1) {
  std::vector<std::string> d;
  for (double day = std::floor(t0 / 86400);
       day <= std::floor(t1 / 86400) && d.size() < 400; ++day)
    d.push_back(utc_date(day * 86400));
  return d;
}

std::string db_day_path(const std::string& db, const std::string& date) {
  return db + "/adsb_" + date + ".plb";
}

bool db_day_bbox(const std::string& plb, double& w, double& s, double& e,
                 double& n) {
  FILE* m = std::fopen((plb + ".meta").c_str(), "r");
  if (!m) return false;
  bool ok = std::fscanf(m, "bbox %lf %lf %lf %lf", &w, &s, &e, &n) == 4;
  std::fclose(m);
  return ok;
}

static std::string shq(const std::string& s) {
  std::string o = "'";
  for (char c : s)
    o += c == '\'' ? std::string("'\\''") : std::string(1, c);
  return o + "'";
}

static bool run_capture(const std::string& cmd, std::string& out) {
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) return false;
  char buf[1 << 14];
  size_t n;
  out.clear();
  while ((n = std::fread(buf, 1, sizeof buf, p)) > 0)
    out.append(buf, n);
  return pclose(p) == 0;
}

static std::string auth_header() {
  const char* t = std::getenv("GITHUB_TOKEN");
  return t && *t ? std::string(" -H ") +
                       shq(std::string("Authorization: Bearer ") + t)
                 : "";
}

// Returns the ordered tar part URLs for a day, trying prod then staging builds.
static std::vector<std::string> release_parts(const std::string& curl,
                                              const std::string& date) {
  std::string tagdate = date;
  for (auto& c : tagdate)
    if (c == '-') c = '.';
  for (const char* flavour : {"prod-0", "prod-1", "staging-0"}) {
    std::string tag = "v" + tagdate + "-planes-readsb-" + flavour;
    std::string url = "https://api.github.com/repos/adsblol/globe_history_" +
                      date.substr(0, 4) + "/releases/tags/" + tag;
    std::string body;
    if (!run_capture(
            curl + " -fsSL" + auth_header() + " " + shq(url) + " 2>/dev/null",
            body))
      continue;
    JVal j;
    if (!parse_json(body, j)) continue;
    const JVal* assets = j.get("assets");
    if (!assets || assets->t != JVal::Arr) continue;
    std::vector<std::pair<std::string, std::string>> parts;
    for (auto& a : assets->a) {
      const JVal *n = a.get("name"), *u = a.get("browser_download_url");
      if (n && u && n->s.find(".tar") != std::string::npos)
        parts.push_back({n->s, u->s});
    }
    if (parts.empty()) continue;  // release created but upload not finished yet
    std::sort(parts.begin(), parts.end());
    std::vector<std::string> urls;
    for (auto& p : parts)
      urls.push_back(p.second);
    log("update %s: release %s (%zu part%s)", date.c_str(), tag.c_str(),
        urls.size(), urls.size() > 1 ? "s" : "");
    return urls;
  }
  return {};
}

std::string latest_available_date(const std::string& curl, int max_back) {
  double now = (double)std::time(nullptr);
  for (int k = 1; k <= max_back; ++k) {
    std::string d = utc_date(now - k * 86400.0);
    if (!release_parts(curl, d).empty()) return d;
  }
  return "";
}

static bool gunzip(const std::string& in, std::string& out) {
  z_stream z{};
  if (inflateInit2(&z, 16 + MAX_WBITS) != Z_OK) return false;
  z.next_in = (Bytef*)in.data();
  z.avail_in = (uInt)in.size();
  out.resize(std::max<size_t>(in.size() * 6, 1 << 16));
  size_t have = 0;
  int rc;
  do {
    if (have == out.size()) out.resize(out.size() * 2);
    z.next_out = (Bytef*)out.data() + have;
    z.avail_out = (uInt)(out.size() - have);
    rc = inflate(&z, Z_NO_FLUSH);
    have = out.size() - z.avail_out;
  } while (rc == Z_OK);
  inflateEnd(&z);
  out.resize(have);
  return rc == Z_STREAM_END;
}

static uint64_t tar_num(const char* p, int n) {
  if ((unsigned char)p[0] & 0x80) {  // GNU base-256
    uint64_t v = 0;
    for (int i = 1; i < n; ++i)
      v = (v << 8) | (unsigned char)p[i];
    return v;
  }
  uint64_t v = 0;
  for (int i = 0; i < n && p[i]; ++i)
    if (p[i] >= '0' && p[i] <= '7') v = v * 8 + (p[i] - '0');
  return v;
}

// Streams `cmd` output as a tar archive; trace JSON members are parsed by a
// worker pool and written to the cache as they arrive.
static bool ingest_tar_stream(const std::string& cmd, const LoadOptions& filt,
                              CacheWriter& w, size_t& nfiles) {
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) return false;
  std::deque<std::string> q;
  std::mutex mu;
  std::condition_variable cv_put, cv_get;
  bool eof = false;
  const size_t kMaxQ = 512;
  int nw = std::max(1, num_threads() - 1);
  std::vector<std::thread> workers;
  std::atomic<size_t> bad{0};
  for (int i = 0; i < nw; ++i)
    workers.emplace_back([&] {
      std::string raw, text;
      std::vector<Flight> fl;
      for (;;) {
        {
          std::unique_lock<std::mutex> lk(mu);
          cv_get.wait(lk, [&] { return !q.empty() || eof; });
          if (q.empty()) return;
          raw = std::move(q.front());
          q.pop_front();
        }
        cv_put.notify_one();
        const std::string* src = &raw;
        if (raw.size() > 2 && (unsigned char)raw[0] == 0x1f &&
            (unsigned char)raw[1] == 0x8b) {
          if (!gunzip(raw, text)) {
            ++bad;
            continue;
          }
          src = &text;
        }
        JVal j;
        if (!parse_json(*src, j)) {
          ++bad;
          continue;
        }
        fl.clear();
        parse_readsb_trace(j, filt, fl);
        for (auto& f : fl) {
          std::sort(
              f.s.begin(), f.s.end(),
              [](const RawSample& a, const RawSample& b) { return a.t < b.t; });
          w.add(f);
        }
      }
    });
  auto readn = [&](char* buf, size_t n) {
    return std::fread(buf, 1, n, p) == n;
  };
  char hdr[512];
  std::string longname;
  bool ok = true;
  Timer tm;
  size_t bytes = 0;
  while (readn(hdr, 512)) {
    bool zero = true;
    for (char c : hdr)
      zero &= c == 0;
    if (zero)
      continue;  // end-of-archive blocks (may repeat between concatenated
                 // parts)
    uint64_t size = tar_num(hdr + 124, 12);
    char type = hdr[156];
    std::string name =
        longname.empty() ? std::string(hdr, strnlen(hdr, 100)) : longname;
    longname.clear();
    if (!std::memcmp(hdr + 257, "ustar", 5) && hdr[345] &&
        name.find('/') == std::string::npos)
      name = std::string(hdr + 345, strnlen(hdr + 345, 155)) + "/" + name;
    uint64_t padded = (size + 511) / 512 * 512;
    std::string data;
    bool want = (type == '0' || type == 0) &&
                name.find("trace_full_") != std::string::npos;
    if (type == 'L' || want) {
      data.resize(padded);
      if (padded && !readn(data.data(), padded)) {
        ok = false;
        break;
      }
      data.resize(size);
    } else {
      std::vector<char> skip(1 << 16);
      for (uint64_t left = padded; left > 0;) {
        size_t k = (size_t)std::min<uint64_t>(left, skip.size());
        if (!readn(skip.data(), k)) {
          ok = false;
          break;
        }
        left -= k;
      }
      if (!ok) break;
    }
    bytes += 512 + padded;
    if (type == 'L') {
      longname.assign(data.c_str());
      continue;
    }
    if (!want) continue;
    ++nfiles;
    {
      std::unique_lock<std::mutex> lk(mu);
      cv_put.wait(lk, [&] { return q.size() < kMaxQ; });
      q.push_back(std::move(data));
    }
    cv_get.notify_one();
    if (g_verbose && nfiles % 5000 == 0)
      std::fprintf(stderr,
                   "  ... %zu traces, %.0f MB, %.0f MB/s, %llu legs kept\r",
                   nfiles, bytes / 1e6, bytes / 1e6 / std::max(1e-3, tm.sec()),
                   (unsigned long long)w.flights());
  }
  {
    std::lock_guard<std::mutex> lk(mu);
    eof = true;
  }
  cv_get.notify_all();
  for (auto& t : workers)
    t.join();
  int rc = pclose(p);
  if (g_verbose) std::fprintf(stderr, "\n");
  if (bad) log("update: %zu unparsable trace files skipped", (size_t)bad);
  return ok && rc == 0;
}

int run_update(const UpdateOptions& o) {
  fs::create_directories(o.db);
  int failures = 0;
  for (const std::string& date : o.dates) {
    std::string out = db_day_path(o.db, date);
    if (!o.force && fs::exists(out)) {
      log("update %s: %s already present (use --force to rebuild)",
          date.c_str(), out.c_str());
      continue;
    }
    auto urls = release_parts(o.curl, date);
    if (urls.empty()) {
      log("update %s: no published adsb.lol release (daily dumps appear ~1 day "
          "later; GitHub API rate limit "
          "is 60/h without GITHUB_TOKEN)",
          date.c_str());
      ++failures;
      continue;
    }
    std::string cmd = "(";
    for (auto& u : urls)
      cmd += o.curl + " -fsSL --retry 3 " + shq(u) + "; ";
    cmd += ")";
    Timer t;
    CacheWriter w(out);
    LoadOptions f = o.filter;
    double d0 = parse_time(date);
    f.tmin = std::max(
        f.tmin, d0 - 3600);  // traces are per UTC day; keep a little overlap
    f.tmax = std::min(f.tmax, d0 + 86400 + 3600);
    size_t nfiles = 0;
    if (!ingest_tar_stream(cmd, f, w, nfiles) || nfiles == 0) {
      log("update %s: download/extract failed after %zu traces", date.c_str(),
          nfiles);
      ++failures;
      continue;  // CacheWriter destructor discards the partial file
    }
    w.close();
    {
      FILE* m = std::fopen((out + ".meta").c_str(), "w");
      if (m) {
        std::fprintf(m,
                     "bbox %.4f %.4f %.4f %.4f\nmax_alt_m %.0f\nsource "
                     "adsb.lol\nbuilt %s\n",
                     f.minlon, f.minlat, f.maxlon, f.maxlat, f.max_alt_m,
                     fmt_time((double)std::time(nullptr)).c_str());
        std::fclose(m);
      }
    }
    log("update %s: %zu traces -> %llu legs, %llu samples -> %s (%.0fs)",
        date.c_str(), nfiles, (unsigned long long)w.flights(),
        (unsigned long long)w.samples(), out.c_str(), t.sec());
  }
  return failures;
}

void print_db(const std::string& db) {
  if (!fs::is_directory(db)) {
    std::printf("database %s: empty (run `planeloc update`)\n", db.c_str());
    return;
  }
  std::vector<std::string> v;
  for (auto& e : fs::directory_iterator(db)) {
    std::string n = e.path().filename().string();
    if (n.rfind("adsb_", 0) == 0 && n.size() > 9 &&
        n.substr(n.size() - 4) == ".plb")
      v.push_back(n.substr(5, 10) + "  " +
                  std::to_string(e.file_size() / (1 << 20)) + " MB");
  }
  std::sort(v.begin(), v.end());
  std::printf("database %s: %zu day(s)\n", db.c_str(), v.size());
  for (auto& s : v)
    std::printf("  %s\n", s.c_str());
}

}  // namespace pl
