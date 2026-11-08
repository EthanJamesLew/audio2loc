// WGS84 geodesy, local frames and the observer search grid.
#pragma once
#include <cmath>
#include <cstdint>

namespace pl {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kWgsA = 6378137.0;
constexpr double kWgsF = 1.0 / 298.257223563;
constexpr double kWgsE2 = kWgsF * (2.0 - kWgsF);
constexpr double kFt = 0.3048;
constexpr double kKt = 0.514444444;

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3() = default;
  Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3& operator+=(const Vec3& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  double norm() const { return std::sqrt(x * x + y * y + z * z); }
  double norm2() const { return x * x + y * y + z * z; }
};

Vec3 geodetic_to_ecef(double lat_deg, double lon_deg, double h);
void ecef_to_geodetic(const Vec3& p, double& lat_deg, double& lon_deg, double& h);
// ENU direction vector at (lat, lon) expressed in ECEF.
Vec3 enu_to_ecef_dir(double lat_deg, double lon_deg, double e, double n, double u);
void ecef_dir_to_enu(double lat_deg, double lon_deg, const Vec3& d, double& e, double& n, double& u);

// Metres per degree of latitude / longitude at a latitude (ellipsoidal).
inline double m_per_deg_lat(double lat_deg) {
  double s = std::sin(lat_deg * kDeg);
  double w = 1.0 - kWgsE2 * s * s;
  return kDeg * kWgsA * (1.0 - kWgsE2) / (w * std::sqrt(w));
}
inline double m_per_deg_lon(double lat_deg) {
  double s = std::sin(lat_deg * kDeg);
  return kDeg * kWgsA * std::cos(lat_deg * kDeg) / std::sqrt(1.0 - kWgsE2 * s * s);
}

double haversine_m(double lat1, double lon1, double lat2, double lon2);

// A quasi-uniform lat/lon grid with ~step_m spacing. Row i sits at latitude i*dlat;
// within a row the longitude spacing is widened by 1/cos(lat) so cells stay ~square.
class GeoGrid {
 public:
  explicit GeoGrid(double step_m);
  double step() const { return step_; }
  static uint64_t pack(int32_t i, int32_t j) { return ((uint64_t)(uint32_t)i << 32) | (uint32_t)j; }
  static int32_t ki(uint64_t k) { return (int32_t)(uint32_t)(k >> 32); }
  static int32_t kj(uint64_t k) { return (int32_t)(uint32_t)(k & 0xffffffffu); }
  double row_lat(int32_t i) const { return i * dlat_; }
  double row_dlon(int32_t i) const;
  void center(uint64_t key, double& lat, double& lon) const;
  uint64_t nearest(double lat, double lon) const;
  // Calls fn(key, lat, lon) for every node within radius_m (horizontal) of (lat, lon).
  template <class F>
  void for_each_within(double lat, double lon, double radius_m, F&& fn) const {
    double mlat = m_per_deg_lat(lat);
    int32_t i0 = (int32_t)std::ceil((lat - radius_m / mlat) / dlat_);
    int32_t i1 = (int32_t)std::floor((lat + radius_m / mlat) / dlat_);
    for (int32_t i = i0; i <= i1; ++i) {
      double la = i * dlat_;
      double dy = (la - lat) * mlat;
      double rem = radius_m * radius_m - dy * dy;
      if (rem < 0) continue;
      double dxm = std::sqrt(rem);
      double mlon = m_per_deg_lon(la);
      double dl = row_dlon(i);
      int32_t j0 = (int32_t)std::ceil((lon - dxm / mlon) / dl);
      int32_t j1 = (int32_t)std::floor((lon + dxm / mlon) / dl);
      for (int32_t j = j0; j <= j1; ++j) fn(pack(i, j), la, j * dl);
    }
  }

 private:
  double step_, dlat_;
};

}  // namespace pl
