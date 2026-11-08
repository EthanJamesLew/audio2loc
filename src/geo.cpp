#include "geo.h"

namespace pl {

Vec3 geodetic_to_ecef(double lat_deg, double lon_deg, double h) {
  double la = lat_deg * kDeg, lo = lon_deg * kDeg;
  double sl = std::sin(la), cl = std::cos(la);
  double N = kWgsA / std::sqrt(1.0 - kWgsE2 * sl * sl);
  return {(N + h) * cl * std::cos(lo), (N + h) * cl * std::sin(lo),
          (N * (1.0 - kWgsE2) + h) * sl};
}

void ecef_to_geodetic(const Vec3& p, double& lat_deg, double& lon_deg,
                      double& h) {
  // Bowring's method with two refinement iterations (sub-mm for aviation
  // altitudes).
  double b = kWgsA * (1.0 - kWgsF);
  double ep2 = (kWgsA * kWgsA - b * b) / (b * b);
  double r = std::sqrt(p.x * p.x + p.y * p.y);
  double lon = std::atan2(p.y, p.x);
  double th = std::atan2(p.z * kWgsA, r * b);
  double lat = std::atan2(p.z + ep2 * b * std::pow(std::sin(th), 3),
                          r - kWgsE2 * kWgsA * std::pow(std::cos(th), 3));
  for (int it = 0; it < 2; ++it) {
    double sl = std::sin(lat);
    double N = kWgsA / std::sqrt(1.0 - kWgsE2 * sl * sl);
    double hh = r / std::cos(lat) - N;
    lat = std::atan2(p.z, r * (1.0 - kWgsE2 * N / (N + hh)));
  }
  double sl = std::sin(lat);
  double N = kWgsA / std::sqrt(1.0 - kWgsE2 * sl * sl);
  h = (std::fabs(std::cos(lat)) > 1e-9) ? r / std::cos(lat) - N
                                        : std::fabs(p.z) - b;
  lat_deg = lat / kDeg;
  lon_deg = lon / kDeg;
}

Vec3 enu_to_ecef_dir(double lat_deg, double lon_deg, double e, double n,
                     double u) {
  double la = lat_deg * kDeg, lo = lon_deg * kDeg;
  double sl = std::sin(la), cl = std::cos(la), so = std::sin(lo),
         co = std::cos(lo);
  return {-so * e - sl * co * n + cl * co * u,
          co * e - sl * so * n + cl * so * u, cl * n + sl * u};
}

void ecef_dir_to_enu(double lat_deg, double lon_deg, const Vec3& d, double& e,
                     double& n, double& u) {
  double la = lat_deg * kDeg, lo = lon_deg * kDeg;
  double sl = std::sin(la), cl = std::cos(la), so = std::sin(lo),
         co = std::cos(lo);
  e = -so * d.x + co * d.y;
  n = -sl * co * d.x - sl * so * d.y + cl * d.z;
  u = cl * co * d.x + cl * so * d.y + sl * d.z;
}

double haversine_m(double lat1, double lon1, double lat2, double lon2) {
  double dla = (lat2 - lat1) * kDeg, dlo = (lon2 - lon1) * kDeg;
  double a = std::sin(dla / 2) * std::sin(dla / 2) +
             std::cos(lat1 * kDeg) * std::cos(lat2 * kDeg) * std::sin(dlo / 2) *
                 std::sin(dlo / 2);
  return 2.0 * 6371008.8 * std::asin(std::sqrt(std::min(1.0, a)));
}

GeoGrid::GeoGrid(double step_m) : step_(step_m), dlat_(step_m / 111132.954) {}

double GeoGrid::row_dlon(int32_t i) const {
  double c = std::cos(i * dlat_ * kDeg);
  if (c < 1e-6) c = 1e-6;
  return step_ / (111319.49 * c);
}

void GeoGrid::center(uint64_t key, double& lat, double& lon) const {
  int32_t i = ki(key), j = kj(key);
  lat = i * dlat_;
  lon = j * row_dlon(i);
}

uint64_t GeoGrid::nearest(double lat, double lon) const {
  int32_t i = (int32_t)std::lround(lat / dlat_);
  int32_t j = (int32_t)std::lround(lon / row_dlon(i));
  return pack(i, j);
}

}  // namespace pl
