#include "dem.h"

#include <cmath>
#include <cstdio>

#include "io.h"
#include "util.h"

namespace pl {

void Dem::preload(double minlat, double maxlat, double minlon, double maxlon) {
  if (dir_.empty()) return;
  int n_loaded = 0;
  for (int la = (int)std::floor(minlat); la <= (int)std::floor(maxlat); ++la) {
    for (int lo = (int)std::floor(minlon); lo <= (int)std::floor(maxlon); ++lo) {
      if (tiles_.count(tkey(la, lo))) continue;
      char name[32];
      std::snprintf(name, sizeof name, "%c%02d%c%03d.hgt", la >= 0 ? 'N' : 'S', std::abs(la), lo >= 0 ? 'E' : 'W',
                    std::abs(lo));
      std::string data;
      bool ok = false;
      for (const char* ext : {"", ".gz"}) {
        std::string p = dir_ + "/" + name + ext;
        if (file_exists(p) && read_file_maybe_gz(p, data)) {
          ok = true;
          break;
        }
      }
      if (!ok) continue;
      Tile t;
      size_t cells = data.size() / 2;
      t.n = (int)std::lround(std::sqrt((double)cells));
      if ((size_t)t.n * t.n != cells) {
        log("dem: %s has unexpected size, skipped", name);
        continue;
      }
      t.h.resize(cells);
      const unsigned char* p = (const unsigned char*)data.data();
      for (size_t i = 0; i < cells; ++i) t.h[i] = (int16_t)((p[2 * i] << 8) | p[2 * i + 1]);  // big endian
      tiles_[tkey(la, lo)] = std::move(t);
      ++n_loaded;
    }
  }
  if (n_loaded) log("dem: loaded %d tile(s) from %s", n_loaded, dir_.c_str());
}

double Dem::elev(double lat, double lon) const {
  if (tiles_.empty()) return def_;
  int la = (int)std::floor(lat), lo = (int)std::floor(lon);
  auto it = tiles_.find(tkey(la, lo));
  if (it == tiles_.end()) return def_;
  const Tile& t = it->second;
  double x = (lon - lo) * (t.n - 1), y = (la + 1 - lat) * (t.n - 1);
  int ix = std::clamp((int)x, 0, t.n - 2), iy = std::clamp((int)y, 0, t.n - 2);
  double fx = x - ix, fy = y - iy;
  auto H = [&](int xx, int yy) {
    int16_t v = t.h[(size_t)yy * t.n + xx];
    return v == -32768 ? def_ : (double)v;  // SRTM void
  };
  return (H(ix, iy) * (1 - fx) + H(ix + 1, iy) * fx) * (1 - fy) + (H(ix, iy + 1) * (1 - fx) + H(ix + 1, iy + 1) * fx) * fy;
}

}  // namespace pl
