// Small shared utilities: logging, threading, time parsing, string helpers.
#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace pl {

extern int g_verbose;
extern int g_threads;

void log(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
void vlog(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
[[noreturn]] void die(const char* fmt, ...)
    __attribute__((format(printf, 1, 2)));

int num_threads();

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// Dynamic-scheduled parallel loop. f(i, thread_id).
template <class F>
void parallel_for(size_t n, F&& f, size_t grain = 1) {
  if (n == 0) return;
  size_t chunks = (n + grain - 1) / grain;
  int nt = (int)std::min<size_t>((size_t)num_threads(), chunks);
  if (nt <= 1) {
    for (size_t i = 0; i < n; ++i)
      f(i, 0);
    return;
  }
  std::atomic<size_t> next{0};
  std::vector<std::thread> th;
  th.reserve(nt);
  for (int t = 0; t < nt; ++t) {
    th.emplace_back([&, t] {
      for (;;) {
        size_t b = next.fetch_add(grain);
        if (b >= n) break;
        size_t e = std::min(n, b + grain);
        for (size_t i = b; i < e; ++i)
          f(i, t);
      }
    });
  }
  for (auto& x : th)
    x.join();
}

// Time: seconds since Unix epoch (UTC), double precision.
double parse_time(const std::string& s);  // epoch number, YYYY-MM-DD,
                                          // YYYY-MM-DDTHH:MM[:SS[.f]][Z|+hh:mm]
std::string fmt_time(double epoch);       // ISO 8601 UTC with milliseconds

std::vector<std::string> split(const std::string& s, char sep);
std::string trim(const std::string& s);
std::string lower(std::string s);
bool parse_double(const std::string& s, double& out);
bool file_exists(const std::string& p);

// Memory accounting (Linux /proc). Bytes.
size_t current_rss_bytes();
size_t mem_available_bytes();
size_t peak_rss_bytes();  // high-water mark of the resident set
// Resolves the memory limit: explicit GB value (>0), else $PLANELOC_MEM_LIMIT
// (GB), else 50% of the memory available right now.
size_t resolve_mem_limit(double gb);
// Background thread that terminates the process cleanly (exit code 3, clear
// message) if the resident set ever exceeds `limit` - a hard backstop so a
// mis-estimate can never take the machine down.
void start_mem_watchdog(size_t limit);
inline double gb(size_t b) {
  return b / 1073741824.0;
}

struct Timer {
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  double sec() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
        .count();
  }
};

inline uint64_t mix64(uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

}  // namespace pl
