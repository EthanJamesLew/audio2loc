#include "adsb.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <unordered_map>

#include "io.h"
#include "util.h"

namespace fs = std::filesystem;

namespace pl {

static bool keep(const LoadOptions& o, double t, double lat, double lon, double alt) {
  if (t < o.tmin || t > o.tmax) return false;
  if (lat < o.minlat || lat > o.maxlat || lon < o.minlon || lon > o.maxlon) return false;
  if (std::isfinite(alt) && alt > o.max_alt_m) return false;
  return true;
}

// ------------------------------------------------------------------ CSV ---
static std::vector<std::string> csv_fields(const std::string& line) {
  std::vector<std::string> out;
  std::string cur;
  bool q = false;
  for (size_t i = 0; i < line.size(); ++i) {
    char c = line[i];
    if (q) {
      if (c == '"' && i + 1 < line.size() && line[i + 1] == '"') cur += '"', ++i;
      else if (c == '"') q = false;
      else cur += c;
    } else if (c == '"') q = true;
    else if (c == ',') out.push_back(cur), cur.clear();
    else if (c != '\r') cur += c;
  }
  out.push_back(cur);
  return out;
}

static void load_csv(const std::string& path, const LoadOptions& o, std::vector<Flight>& out) {
  std::string text;
  if (!read_file_maybe_gz(path, text)) die("cannot read %s", path.c_str());
  size_t pos = 0;
  auto next_line = [&](std::string& l) {
    if (pos >= text.size()) return false;
    size_t e = text.find('\n', pos);
    if (e == std::string::npos) e = text.size();
    l.assign(text, pos, e - pos);
    pos = e + 1;
    return true;
  };
  std::string line;
  if (!next_line(line)) return;
  auto hdr = csv_fields(line);
  std::unordered_map<std::string, int> col;
  for (size_t i = 0; i < hdr.size(); ++i) col[lower(trim(hdr[i]))] = (int)i;
  auto find = [&](std::initializer_list<const char*> names) {
    for (auto n : names) {
      auto it = col.find(n);
      if (it != col.end()) return it->second;
    }
    return -1;
  };
  int c_t = find({"time", "timestamp", "ts", "epoch", "datetime"});
  int c_id = find({"icao24", "icao", "hex", "icao_address", "addr"});
  int c_cs = find({"callsign", "flight", "ident"});
  int c_ty = find({"typecode", "type", "aircraft_type", "actype", "t"});
  int c_la = find({"lat", "latitude"});
  int c_lo = find({"lon", "lng", "longitude"});
  int c_gnd = find({"onground", "on_ground", "ground"});
  // Altitude: prefer geometric. Unit detection: explicit suffixes, OpenSky raw = metres, traffic lib = feet.
  struct ColU { int c; double scale; bool geom; };
  std::vector<std::pair<const char*, ColU>> alt_cands = {
      {"alt_geom_m", {-1, 1, true}},   {"alt_m", {-1, 1, false}},       {"alt_geom_ft", {-1, kFt, true}},
      {"alt_geom", {-1, kFt, true}},   {"geoaltitude", {-1, 1, true}},  {"alt_ft", {-1, kFt, false}},
      {"alt_baro", {-1, kFt, false}},  {"baroaltitude", {-1, 1, false}}, {"altitude", {-1, kFt, false}},
      {"alt", {-1, kFt, false}}};
  ColU alt{-1, 1, false};
  for (auto& [n, cu] : alt_cands) {
    auto it = col.find(n);
    if (it != col.end()) {
      alt = cu;
      alt.c = it->second;
      break;
    }
  }
  if (alt.c >= 0 && !o.csv_alt_unit.empty()) alt.scale = o.csv_alt_unit == "ft" ? kFt : 1.0;
  ColU gs{-1, 1, false};
  for (auto [n, s] : std::vector<std::pair<const char*, double>>{
           {"gs_mps", 1}, {"velocity", 1}, {"gs_kt", kKt}, {"gs", kKt}, {"groundspeed", kKt}, {"speed", kKt}}) {
    auto it = col.find(n);
    if (it != col.end()) {
      gs = {it->second, s, false};
      break;
    }
  }
  if (gs.c >= 0 && !o.csv_speed_unit.empty()) gs.scale = o.csv_speed_unit == "kt" ? kKt : 1.0;
  int c_trk = find({"track", "heading", "true_track", "trk"});
  ColU vr{-1, 1, false};
  for (auto [n, s] : std::vector<std::pair<const char*, double>>{{"vrate_mps", 1},
                                                                 {"vertrate", 1},
                                                                 {"vrate_fpm", kFt / 60},
                                                                 {"vertical_rate", kFt / 60},
                                                                 {"baro_rate", kFt / 60},
                                                                 {"geom_rate", kFt / 60}}) {
    auto it = col.find(n);
    if (it != col.end()) {
      vr = {it->second, s, false};
      break;
    }
  }
  if (c_t < 0 || c_id < 0 || c_la < 0 || c_lo < 0 || alt.c < 0)
    die("%s: CSV needs columns time, icao24, lat, lon and an altitude column (alt_m/alt_ft/geoaltitude/...)",
        path.c_str());
  std::unordered_map<std::string, size_t> idx;
  size_t base = out.size();
  long nrows = 0, nkept = 0;
  auto num = [](const std::string& s) {
    double v;
    return parse_double(s, v) ? v : kNaN;
  };
  while (next_line(line)) {
    if (line.empty()) continue;
    auto f = csv_fields(line);
    ++nrows;
    if ((int)f.size() < (int)hdr.size()) continue;
    if (c_gnd >= 0) {
      std::string g = lower(trim(f[c_gnd]));
      if (g == "true" || g == "1" || g == "t") continue;
    }
    double t;
    try {
      t = parse_time(f[c_t]);
    } catch (...) {
      continue;
    }
    double la = num(f[c_la]), lo = num(f[c_lo]), al = num(f[alt.c]) * alt.scale;
    if (!std::isfinite(la) || !std::isfinite(lo) || !std::isfinite(al)) continue;
    if (alt.geom) al -= o.geoid_undulation;
    if (!keep(o, t, la, lo, al)) continue;
    std::string id = lower(trim(f[c_id]));
    auto it = idx.find(id);
    if (it == idx.end()) {
      it = idx.emplace(id, out.size()).first;
      out.emplace_back();
      out.back().icao = id;
    }
    Flight& fl = out[it->second];
    if (c_cs >= 0 && fl.callsign.empty()) fl.callsign = trim(f[c_cs]);
    if (c_ty >= 0 && fl.type.empty()) fl.type = trim(f[c_ty]);
    RawSample s{t, la, lo, (float)al, gs.c >= 0 ? (float)(num(f[gs.c]) * gs.scale) : NAN,
                c_trk >= 0 ? (float)num(f[c_trk]) : NAN, vr.c >= 0 ? (float)(num(f[vr.c]) * vr.scale) : NAN};
    fl.s.push_back(s);
    ++nkept;
  }
  for (size_t i = base; i < out.size(); ++i) out[i].cls = classify_type(out[i].type);
  vlog("csv %s: %ld rows, %ld kept, %zu aircraft", path.c_str(), nrows, nkept, out.size() - base);
}

// --------------------------------------------- readsb / adsb.lol traces ---
// Format: {"icao":..,"r":..,"t":type,"timestamp":base,"trace":[[dt,lat,lon,alt_ft|"ground",gs_kt,track,flags,
//          vrate_fpm,{details}|null,source,geom_alt_ft,geom_rate,ias,roll],...]}
void parse_readsb_trace(const JVal& root, const LoadOptions& o, std::vector<Flight>& out) {
  const JVal* tr = root.get("trace");
  const JVal* ts = root.get("timestamp");
  if (!tr || tr->t != JVal::Arr || !ts || !ts->is_num()) return;
  Flight cur;
  if (auto* v = root.get("icao")) cur.icao = lower(v->s);
  if (auto* v = root.get("t"); v && v->t == JVal::Str) cur.type = v->s;
  cur.cls = classify_type(cur.type);
  auto flush = [&] {
    if (cur.s.size() >= 2) out.push_back(cur);
    cur.s.clear();
  };
  for (const JVal& p : tr->a) {
    if (p.t != JVal::Arr || p.a.size() < 6) continue;
    int flags = p.a.size() > 6 && p.a[6].is_num() ? (int)p.a[6].n : 0;
    if (flags & 2) flush();  // new leg
    if (p.a.size() > 8 && p.a[8].t == JVal::Obj && cur.callsign.empty())
      if (auto* fl = p.a[8].get("flight")) cur.callsign = trim(fl->s);
    if (!p.a[3].is_num()) continue;  // "ground" or null altitude
    double t = ts->n + p.a[0].n, la = p.a[1].n, lo = p.a[2].n;
    double alt = p.a[3].n * kFt;
    if (p.a.size() > 10 && p.a[10].is_num()) alt = p.a[10].n * kFt - o.geoid_undulation;
    if (!keep(o, t, la, lo, alt)) continue;
    RawSample s{t,
                la,
                lo,
                (float)alt,
                p.a[4].is_num() ? (float)(p.a[4].n * kKt) : NAN,
                p.a[5].is_num() ? (float)p.a[5].n : NAN,
                p.a.size() > 7 && p.a[7].is_num() ? (float)(p.a[7].n * kFt / 60) : NAN};
    cur.s.push_back(s);
  }
  flush();
}

// ---------------------------------------------------------------- cache ---
static const char kMagic[8] = {'P', 'L', 'B', 'C', 'A', 'C', 'H', '1'};

static void write_flight(FILE* f, const Flight& x) {
  auto ws = [&](const std::string& s) {
    uint32_t l = (uint32_t)s.size();
    std::fwrite(&l, 4, 1, f);
    std::fwrite(s.data(), 1, l, f);
  };
  ws(x.icao);
  ws(x.callsign);
  ws(x.type);
  uint64_t m = x.s.size();
  std::fwrite(&m, 8, 1, f);
  std::fwrite(x.s.data(), sizeof(RawSample), m, f);
}

void save_cache(const std::string& path, const std::vector<Flight>& fl) {
  CacheWriter w(path);
  for (auto& x : fl) w.add(x);
  w.close();
}

CacheWriter::CacheWriter(const std::string& path) : path_(path), tmp_(path + ".part") {
  f_ = std::fopen(tmp_.c_str(), "wb");
  if (!f_) die("cannot write %s", tmp_.c_str());
  std::fwrite(kMagic, 1, 8, f_);
  uint64_t n = 0;
  std::fwrite(&n, 8, 1, f_);  // patched in close()
}

void CacheWriter::add(const Flight& f) {
  std::lock_guard<std::mutex> g(mu_);
  write_flight(f_, f);
  ++n_;
  samples_ += f.s.size();
}

void CacheWriter::close() {
  if (!f_) return;
  std::fseek(f_, 8, SEEK_SET);
  std::fwrite(&n_, 8, 1, f_);
  bool ok = std::fflush(f_) == 0;
  std::fclose(f_);
  f_ = nullptr;
  if (!ok || std::rename(tmp_.c_str(), path_.c_str()) != 0) die("failed to finalise %s", path_.c_str());
}

CacheWriter::~CacheWriter() {
  if (f_) {  // not closed: abandon the partial file
    std::fclose(f_);
    std::remove(tmp_.c_str());
  }
}

static void load_cache(const std::string& path, const LoadOptions& o, std::vector<Flight>& out) {
  FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) die("cannot read %s", path.c_str());
  char mg[8];
  uint64_t n = 0;
  if (std::fread(mg, 1, 8, f) != 8 || std::memcmp(mg, kMagic, 8) || std::fread(&n, 8, 1, f) != 1)
    die("%s: not a planeloc cache", path.c_str());
  auto rs = [&](std::string& s) {
    uint32_t l = 0;
    if (std::fread(&l, 4, 1, f) != 1) die("truncated cache");
    s.resize(l);
    if (l && std::fread(s.data(), 1, l, f) != l) die("truncated cache");
  };
  for (uint64_t i = 0; i < n; ++i) {
    Flight x;
    rs(x.icao);
    rs(x.callsign);
    rs(x.type);
    uint64_t m = 0;
    if (std::fread(&m, 8, 1, f) != 1) die("truncated cache");
    std::vector<RawSample> s(m);
    if (m && std::fread(s.data(), sizeof(RawSample), m, f) != m) die("truncated cache");
    double last = -1e300;
    for (auto& r : s)
      if (r.t - last >= o.min_dt && keep(o, r.t, r.lat, r.lon, r.alt_m)) {
        x.s.push_back(r);
        last = r.t;
      }
    x.s.shrink_to_fit();
    x.cls = classify_type(x.type);
    if (x.s.size() >= 2) out.push_back(std::move(x));
  }
  std::fclose(f);
}

static bool ends_with(const std::string& s, const char* suf) {
  size_t n = std::strlen(suf);
  return s.size() >= n && s.compare(s.size() - n, n, suf) == 0;
}

void load_adsb(const std::string& path, const LoadOptions& o, std::vector<Flight>& out) {
  Timer tm;
  size_t before = out.size();
  if (fs::is_directory(path)) {
    std::vector<std::string> files;
    for (auto& e : fs::recursive_directory_iterator(path)) {
      if (!e.is_regular_file()) continue;
      std::string p = e.path().string();
      if (ends_with(p, ".json") || ends_with(p, ".json.gz") || ends_with(p, ".csv") || ends_with(p, ".csv.gz") ||
          ends_with(p, ".plb"))
        files.push_back(p);
    }
    log("adsb: scanning %zu files in %s", files.size(), path.c_str());
    int nt = num_threads();
    std::vector<std::vector<Flight>> part(nt);
    std::mutex mu;
    parallel_for(files.size(), [&](size_t i, int tid) {
      const std::string& p = files[i];
      if (ends_with(p, ".csv") || ends_with(p, ".csv.gz") || ends_with(p, ".plb")) {
        std::vector<Flight> tmp;
        if (ends_with(p, ".plb")) load_cache(p, o, tmp);
        else load_csv(p, o, tmp);
        for (auto& x : tmp) part[tid].push_back(std::move(x));
        return;
      }
      std::string text;
      if (!read_file_maybe_gz(p, text)) return;
      JVal root;
      if (!parse_json(text, root)) {
        vlog("adsb: skipping unparsable %s", p.c_str());
        return;
      }
      parse_readsb_trace(root, o, part[tid]);
    });
    for (auto& v : part)
      for (auto& x : v) out.push_back(std::move(x));
  } else if (ends_with(path, ".plb")) {
    load_cache(path, o, out);
  } else if (ends_with(path, ".json") || ends_with(path, ".json.gz")) {
    std::string text;
    JVal root;
    if (!read_file_maybe_gz(path, text) || !parse_json(text, root)) die("cannot parse %s", path.c_str());
    parse_readsb_trace(root, o, out);
  } else {
    load_csv(path, o, out);
  }
  for (size_t i = before; i < out.size(); ++i) {
    auto& s = out[i].s;
    std::sort(s.begin(), s.end(), [](const RawSample& a, const RawSample& b) { return a.t < b.t; });
    if (o.min_dt > 0 && s.size() > 1) {
      size_t w = 1;
      for (size_t r = 1; r < s.size(); ++r)
        if (s[r].t - s[w - 1].t >= o.min_dt) s[w++] = s[r];
      s.resize(w);
      s.shrink_to_fit();
    }
  }
  size_t ns = 0;
  for (size_t i = before; i < out.size(); ++i) ns += out[i].s.size();
  log("adsb: %s -> %zu flight legs, %zu samples (%.1fs)", path.c_str(), out.size() - before, ns, tm.sec());
}

// --------------------------------------------------------------- tracks ---
void Track::eval(double t, Vec3& p, Vec3& v, int* hint) const {
  int n = (int)pts.size();
  if (t <= pts[0].t) {
    p = pts[0].p + pts[0].v * (t - pts[0].t);
    v = pts[0].v;
    return;
  }
  if (t >= pts[n - 1].t) {
    p = pts[n - 1].p + pts[n - 1].v * (t - pts[n - 1].t);
    v = pts[n - 1].v;
    return;
  }
  int i;
  if (hint && *hint >= 0 && *hint < n - 1) {
    i = *hint;
    while (i < n - 2 && pts[i + 1].t <= t) ++i;
    while (i > 0 && pts[i].t > t) --i;
  } else {
    i = (int)(std::upper_bound(pts.begin(), pts.end(), t, [](double x, const TrackPoint& q) { return x < q.t; }) -
              pts.begin()) - 1;
  }
  if (hint) *hint = i;
  const TrackPoint &a = pts[i], &b = pts[i + 1];
  double h = b.t - a.t, s = (t - a.t) / h, s2 = s * s, s3 = s2 * s;
  double h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s, h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
  double d00 = 6 * s2 - 6 * s, d10 = 3 * s2 - 4 * s + 1, d01 = -6 * s2 + 6 * s, d11 = 3 * s2 - 2 * s;
  p = a.p * h00 + a.v * (h10 * h) + b.p * h01 + b.v * (h11 * h);
  v = (a.p * d00 + a.v * (d10 * h) + b.p * d01 + b.v * (d11 * h)) * (1.0 / h);
}

std::vector<Track> build_tracks(const std::vector<Flight>& flights, double max_gap) {
  const double kMinDt = 3.0;  // decimate dense reports: velocity comes from reported gs/track, not jitter
  int nt = num_threads();
  std::vector<std::vector<Track>> part(nt);
  parallel_for(flights.size(), [&](size_t fi, int tid) {
    const Flight& F = flights[fi];
    // 1) decimate + drop outliers (implied speed > 400 m/s)
    std::vector<const RawSample*> sel;
    for (const auto& s : F.s) {
      if (!sel.empty()) {
        const RawSample* p = sel.back();
        double dt = s.t - p->t;
        if (dt < kMinDt) continue;
        double d = haversine_m(p->lat, p->lon, s.lat, s.lon);
        if (dt < max_gap && d / dt > 400) continue;
      }
      sel.push_back(&s);
    }
    // 2) split at gaps
    size_t a = 0;
    while (a < sel.size()) {
      size_t b = a + 1;
      while (b < sel.size() && sel[b]->t - sel[b - 1]->t <= max_gap) ++b;
      if (b - a >= 2) {
        Track T;
        T.flight = (uint32_t)fi;
        T.cls = F.cls;
        size_t n = b - a;
        T.pts.resize(n);
        for (size_t k = 0; k < n; ++k) {
          const RawSample* s = sel[a + k];
          T.pts[k].t = s->t;
          T.pts[k].p = geodetic_to_ecef(s->lat, s->lon, s->alt_m);
          T.pts[k].lat = (float)s->lat;
          T.pts[k].lon = (float)s->lon;
          T.minlat = std::min(T.minlat, s->lat);
          T.maxlat = std::max(T.maxlat, s->lat);
          T.minlon = std::min(T.minlon, s->lon);
          T.maxlon = std::max(T.maxlon, s->lon);
        }
        for (size_t k = 0; k < n; ++k) {
          const RawSample* s = sel[a + k];
          size_t k0 = k ? k - 1 : k, k1 = k + 1 < n ? k + 1 : k;
          Vec3 fd = (T.pts[k1].p - T.pts[k0].p) * (1.0 / (T.pts[k1].t - T.pts[k0].t));
          double fe, fn, fu;
          ecef_dir_to_enu(s->lat, s->lon, fd, fe, fn, fu);
          double e = fe, nn = fn, u = fu;
          if (std::isfinite(s->gs_mps) && std::isfinite(s->track_deg)) {
            double re = s->gs_mps * std::sin(s->track_deg * kDeg), rn = s->gs_mps * std::cos(s->track_deg * kDeg);
            double span = T.pts[k1].t - T.pts[k0].t;
            if (std::hypot(re - fe, rn - fn) < 30 || span > 30) e = re, nn = rn;
          }
          if (std::isfinite(s->vrate_mps)) u = s->vrate_mps;
          T.pts[k].v = enu_to_ecef_dir(s->lat, s->lon, e, nn, u);
        }
        T.t0 = T.pts.front().t;
        T.t1 = T.pts.back().t;
        part[tid].push_back(std::move(T));
      }
      a = b;
    }
  });
  std::vector<Track> out;
  for (auto& v : part)
    for (auto& t : v) out.push_back(std::move(t));
  std::sort(out.begin(), out.end(), [](const Track& x, const Track& y) {
    return x.flight != y.flight ? x.flight < y.flight : x.t0 < y.t0;
  });
  return out;
}

std::string track_label(const std::vector<Flight>& flights, const Track& tr) {
  const Flight& f = flights[tr.flight];
  std::string s = f.icao;
  if (!f.callsign.empty()) s += "/" + f.callsign;
  if (!f.type.empty()) s += "/" + f.type;
  return s;
}

}  // namespace pl
