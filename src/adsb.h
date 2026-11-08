// ADS-B ingestion and continuous 3D trajectories.
#pragma once
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <string>
#include <vector>

#include "geo.h"

namespace pl {

enum class ACClass : uint8_t { Unknown = 0, Jet = 1, Turboprop = 2, Piston = 3, Heli = 4 };
const char* class_name(ACClass c);
ACClass classify_type(const std::string& icao_type);
// Optional user table: CSV "typecode,class" where class in {jet,turboprop,piston,heli}.
void load_type_table(const std::string& path);

struct RawSample {
  double t;           // epoch s
  double lat, lon;    // deg
  float alt_m;        // metres above MSL (NaN if unknown)
  float gs_mps;       // ground speed (NaN if unknown)
  float track_deg;    // true track (NaN if unknown)
  float vrate_mps;    // vertical rate (NaN if unknown)
};

struct Flight {
  std::string icao, callsign, type;
  ACClass cls = ACClass::Unknown;
  std::vector<RawSample> s;  // time-sorted, airborne only
};

struct LoadOptions {
  double tmin = -1e300, tmax = 1e300;
  double minlat = -90, maxlat = 90, minlon = -180, maxlon = 180;
  double max_alt_m = 7000;        // drop samples above this (inaudible cruise traffic)
  double geoid_undulation = 0;    // N (m): MSL = HAE - N, applied to geometric altitudes
  std::string csv_alt_unit;       // "", "m", "ft" (override auto detection)
  std::string csv_speed_unit;     // "", "mps", "kt"
  double min_dt = 3.0;            // thin samples closer than this (s) at load time; = track decimation step
};

// Loads a path: a .csv/.csv.gz file, a .plb cache, a readsb/adsb.lol trace JSON (gz ok),
// or a directory (recursively: trace_full_*.json, *.json, *.json.gz, *.csv).
void load_adsb(const std::string& path, const LoadOptions& o, std::vector<Flight>& out);
void save_cache(const std::string& path, const std::vector<Flight>& f);

struct JVal;
// Parses one readsb / adsb.lol trace_full JSON document (new-leg flags split it into legs).
void parse_readsb_trace(const JVal& root, const LoadOptions& o, std::vector<Flight>& out);

// Streaming, thread-safe .plb writer (constant memory). Written to <path>.part, renamed on close().
class CacheWriter {
 public:
  explicit CacheWriter(const std::string& path);
  ~CacheWriter();
  void add(const Flight& f);
  void close();
  uint64_t flights() const { return n_; }
  uint64_t samples() const { return samples_; }

 private:
  std::string path_, tmp_;
  FILE* f_ = nullptr;
  std::mutex mu_;
  uint64_t n_ = 0, samples_ = 0;
};

// ---------------------------------------------------------------------------
struct TrackPoint {
  double t;
  Vec3 p, v;  // ECEF position (height = MSL) and velocity
  float lat, lon;  // geodetic position (deg), for fast spatial culling
};

// A contiguous airborne leg with cubic-Hermite interpolation of position, whose derivative
// gives a velocity consistent with the path (heading changes, accelerations, climbs).
struct Track {
  uint32_t flight = 0;
  ACClass cls = ACClass::Unknown;
  double t0 = 0, t1 = 0;
  double minlat = 90, maxlat = -90, minlon = 180, maxlon = -180;
  std::vector<TrackPoint> pts;
  // Evaluate at time t (clamped to span). `hint` speeds up sequential access.
  void eval(double t, Vec3& p, Vec3& v, int* hint = nullptr) const;
  bool covers(double t) const { return t >= t0 && t <= t1; }
};

std::vector<Track> build_tracks(const std::vector<Flight>& flights, double max_gap_s);
std::string track_label(const std::vector<Flight>& flights, const Track& tr);

}  // namespace pl
