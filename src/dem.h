// Digital elevation model: SRTM/Copernicus-style .hgt tiles (1" or 3"),
// bilinear sampling. Elevations are metres above mean sea level. Missing tiles
// fall back to a constant.
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace pl {

class Dem {
 public:
  Dem(std::string dir, double default_elev)
      : dir_(std::move(dir)), def_(default_elev) {}
  // Load every tile intersecting the bounding box. After this, elev() is
  // lock-free / thread-safe.
  void preload(double minlat, double maxlat, double minlon, double maxlon);
  double elev(double lat, double lon) const;
  bool has_data() const { return !tiles_.empty(); }
  double default_elev() const { return def_; }

  struct Tile {
    int n = 0;
    std::vector<int16_t> h;  // n*n, row 0 = north edge
  };
  // key = lat_floor * 1000 + lon_floor (see tkey); exposed for GPU upload
  const std::map<int64_t, Tile>& tiles() const { return tiles_; }

 private:
  std::string dir_;
  double def_;
  std::map<int64_t, Tile> tiles_;
  static int64_t tkey(int la, int lo) { return (int64_t)la * 1000 + lo; }
};

}  // namespace pl
