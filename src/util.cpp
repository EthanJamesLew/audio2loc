#include "util.h"

#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#include <mach/mach.h>
#include <sys/sysctl.h>
#endif

#include <filesystem>
#include <fstream>

#include <cctype>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ctime>

namespace pl {

int g_verbose = 1;
int g_threads = 0;

static const auto g_t0 = std::chrono::steady_clock::now();
static void stamp() {
  static const bool on = [] {
    const char* e = std::getenv("PLANELOC_TIMESTAMPS");
    return e && *e && *e != '0';
  }();
  if (on)
    std::fprintf(
        stderr, "[%8.3f] ",
        std::chrono::duration<double>(std::chrono::steady_clock::now() - g_t0)
            .count());
}

void log(const char* fmt, ...) {
  if (g_verbose < 1) return;
  stamp();
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
  std::fputc('\n', stderr);
}

void vlog(const char* fmt, ...) {
  if (g_verbose < 2) return;
  stamp();
  va_list ap;
  va_start(ap, fmt);
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
  std::fputc('\n', stderr);
}

void die(const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  std::fprintf(stderr, "error: ");
  std::vfprintf(stderr, fmt, ap);
  va_end(ap);
  std::fputc('\n', stderr);
  std::exit(1);
}

int num_threads() {
  if (g_threads > 0) return g_threads;
  unsigned n = std::thread::hardware_concurrency();
  return n ? (int)n : 4;
}

bool parse_double(const std::string& s, double& out) {
  std::string t = trim(s);
  if (t.empty()) return false;
  char* end = nullptr;
  out = std::strtod(t.c_str(), &end);
  return end && *end == '\0';
}

double parse_time(const std::string& in) {
  std::string s = trim(in);
  double num;
  if (parse_double(s, num)) return num;
  int Y = 0, M = 0, D = 0, h = 0, m = 0;
  double sec = 0;
  if (std::sscanf(s.c_str(), "%d-%d-%d", &Y, &M, &D) != 3)
    throw std::runtime_error("bad time: " + s);
  size_t p = 10;
  if (p < s.size() && (s[p] == 'T' || s[p] == ' ' || s[p] == 't')) {
    ++p;
    int n = 0;
    if (std::sscanf(s.c_str() + p, "%d:%d%n", &h, &m, &n) < 2)
      throw std::runtime_error("bad time: " + s);
    p += n;
    if (p < s.size() && s[p] == ':') {
      ++p;
      char* end;
      sec = std::strtod(s.c_str() + p, &end);
      p = end - s.c_str();
    }
  }
  double tz = 0;
  if (p < s.size()) {
    char c = s[p];
    if (c == 'Z' || c == 'z') {
    } else if (c == '+' || c == '-') {
      int th = 0, tm = 0;
      std::sscanf(s.c_str() + p + 1, "%d:%d", &th, &tm);
      if (s.size() - p - 1 == 4 &&
          s.find(':', p) == std::string::npos) {  // +hhmm
        int v = std::atoi(s.c_str() + p + 1);
        th = v / 100;
        tm = v % 100;
      }
      tz = (th * 3600 + tm * 60) * (c == '+' ? 1 : -1);
    }
  }
  std::tm tmv{};
  tmv.tm_year = Y - 1900;
  tmv.tm_mon = M - 1;
  tmv.tm_mday = D;
  tmv.tm_hour = h;
  tmv.tm_min = m;
  tmv.tm_sec = 0;
  return (double)timegm(&tmv) + sec - tz;
}

std::string fmt_time(double t) {
  if (!std::isfinite(t)) return "nan";
  time_t s = (time_t)std::floor(t);
  int ms = (int)std::lround((t - (double)s) * 1000.0);
  if (ms >= 1000) {
    ms -= 1000;
    ++s;
  }
  std::tm tmv{};
  gmtime_r(&s, &tmv);
  char buf[64];
  std::snprintf(buf, sizeof buf, "%04d-%02d-%02dT%02d:%02d:%02d.%03dZ",
                tmv.tm_year + 1900, tmv.tm_mon + 1, tmv.tm_mday, tmv.tm_hour,
                tmv.tm_min, tmv.tm_sec, ms);
  return buf;
}

std::vector<std::string> split(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == sep) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  out.push_back(cur);
  return out;
}

std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace((unsigned char)s[a]))
    ++a;
  while (b > a && std::isspace((unsigned char)s[b - 1]))
    --b;
  return s.substr(a, b - a);
}

std::string lower(std::string s) {
  for (auto& c : s)
    c = (char)std::tolower((unsigned char)c);
  return s;
}

#if defined(__APPLE__)
size_t current_rss_bytes() {
  mach_task_basic_info info{};
  mach_msg_type_number_t n = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&info), &n) != KERN_SUCCESS)
    return 0;
  return info.resident_size;
}

size_t peak_rss_bytes() {
  rusage u{};
  getrusage(RUSAGE_SELF, &u);
  return static_cast<size_t>(u.ru_maxrss);  // bytes on macOS
}

// macOS has no MemAvailable; free plus inactive pages is the usual stand-in.
size_t mem_available_bytes() {
  vm_statistics64_data_t vm{};
  mach_msg_type_number_t n = HOST_VM_INFO64_COUNT;
  if (host_statistics64(mach_host_self(), HOST_VM_INFO64,
                        reinterpret_cast<host_info64_t>(&vm),
                        &n) != KERN_SUCCESS)
    return size_t{8} << 30;
  return (static_cast<size_t>(vm.free_count) + vm.inactive_count) *
         static_cast<size_t>(sysconf(_SC_PAGESIZE));
}

std::string self_exe_path() {
  uint32_t n = 0;
  _NSGetExecutablePath(nullptr, &n);
  std::string p(n, '\0');
  if (_NSGetExecutablePath(p.data(), &n) != 0) return "planeloc";
  p.resize(std::strlen(p.c_str()));
  std::error_code ec;
  auto c = std::filesystem::canonical(p, ec);
  return ec ? p : c.string();
}
#else
size_t current_rss_bytes() {
  long pages = 0, rss = 0;
  FILE* f = std::fopen("/proc/self/statm", "r");
  if (!f) return 0;
  if (std::fscanf(f, "%ld %ld", &pages, &rss) != 2) rss = 0;
  std::fclose(f);
  return (size_t)rss * (size_t)sysconf(_SC_PAGESIZE);
}

size_t peak_rss_bytes() {
  std::ifstream f("/proc/self/status");
  std::string k;
  while (f >> k) {
    if (k == "VmHWM:") {
      size_t v;
      f >> v;
      return v * 1024;
    }
  }
  return 0;
}

size_t mem_available_bytes() {
  std::ifstream f("/proc/meminfo");
  std::string k;
  size_t v;
  std::string unit;
  while (f >> k >> v >> unit)
    if (k == "MemAvailable:") return v * 1024;
  return (size_t)8 << 30;
}

std::string self_exe_path() {
  std::error_code ec;
  auto p = std::filesystem::read_symlink("/proc/self/exe", ec);
  return ec ? std::string("planeloc") : p.string();
}
#endif

// Tightest cgroup-v2 memory.max on the path of this process (containers,
// systemd scopes), or 0.
static size_t cgroup_limit_bytes() {
  std::ifstream f("/proc/self/cgroup");
  std::string line, path;
  while (std::getline(f, line))
    if (line.rfind("0::", 0) == 0) path = line.substr(3);
  size_t best = 0;
  for (;;) {
    std::ifstream m("/sys/fs/cgroup" + path + "/memory.max");
    std::string v;
    if (m >> v && v != "max") {
      size_t b = std::strtoull(v.c_str(), nullptr, 10);
      if (b && (!best || b < best)) best = b;
    }
    if (path.empty() || path == "/") break;
    size_t p = path.find_last_of('/');
    path =
        p == std::string::npos || p == 0 ? std::string("/") : path.substr(0, p);
  }
  return best;
}

size_t resolve_mem_limit(double g) {
  if (g <= 0) {
    const char* e = std::getenv("PLANELOC_MEM_LIMIT");
    if (e && *e) g = std::atof(e);
  }
  size_t cg = cgroup_limit_bytes();
  size_t lim = g > 0 ? (size_t)(g * 1073741824.0)
                     : current_rss_bytes() + mem_available_bytes() / 2;
  if (cg && lim > cg * 85 / 100)
    lim = cg * 85 / 100;  // stay under a container / scope limit
  return lim;
}

void start_mem_watchdog(size_t limit) {
  std::thread([limit] {
    for (;;) {
      size_t r = current_rss_bytes();
      if (r > limit) {
        std::fprintf(stderr,
                     "\nerror: memory limit exceeded (%.1f GB resident > %.1f "
                     "GB limit). Narrow the search "
                     "(--bbox, shorter --start/--end, --rec-start) or raise "
                     "--mem-limit.\n",
                     gb(r), gb(limit));
        std::fflush(stderr);
        _exit(3);
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(150));
    }
  }).detach();
}

bool file_exists(const std::string& p) {
  struct stat st;
  return ::stat(p.c_str(), &st) == 0;
}

}  // namespace pl
