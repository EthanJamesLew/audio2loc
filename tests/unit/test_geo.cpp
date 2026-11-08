#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "geo.h"

using Catch::Approx;
using namespace pl;

TEST_CASE("geodetic and ECEF round-trip", "[geo]") {
  for (double lat : {-89.0, -33.9, 0.0, 38.55, 89.0}) {
    for (double lon : {-179.0, -121.9, 0.0, 151.2}) {
      const Vec3 p = geodetic_to_ecef(lat, lon, 1234.5);
      double la = 0, lo = 0, h = 0;
      ecef_to_geodetic(p, la, lo, h);
      CHECK(la == Approx(lat).margin(1e-9));
      CHECK(lo == Approx(lon).margin(1e-9));
      CHECK(h == Approx(1234.5).margin(1e-4));
    }
  }
}

TEST_CASE("ENU directions survive a round trip through ECEF", "[geo]") {
  const Vec3 d = enu_to_ecef_dir(40.8, -73.4, 3.0, -4.0, 12.0);
  double e = 0, n = 0, u = 0;
  ecef_dir_to_enu(40.8, -73.4, d, e, n, u);
  CHECK(e == Approx(3.0));
  CHECK(n == Approx(-4.0));
  CHECK(u == Approx(12.0));
  CHECK(d.norm() == Approx(13.0));
}

TEST_CASE("haversine matches the metres-per-degree scale", "[geo]") {
  // One degree of latitude is about 111 km; the spherical and ellipsoidal
  // figures agree to well under one percent.
  const double d = haversine_m(40.0, -73.0, 41.0, -73.0);
  CHECK(d == Approx(m_per_deg_lat(40.5)).epsilon(0.005));
  CHECK(haversine_m(12.3, 45.6, 12.3, 45.6) == Approx(0.0).margin(1e-9));
}

TEST_CASE("GeoGrid nodes are about one step apart", "[geo]") {
  const GeoGrid g(100.0);
  const uint64_t k = g.nearest(38.5466, -121.9270);
  double lat = 0, lon = 0;
  g.center(k, lat, lon);
  // The nearest node lies within half a cell diagonal.
  CHECK(haversine_m(lat, lon, 38.5466, -121.9270) < 100.0 * 0.75);
  CHECK(GeoGrid::ki(k) == GeoGrid::ki(GeoGrid::pack(GeoGrid::ki(k), -7)));
  CHECK(GeoGrid::kj(GeoGrid::pack(3, -7)) == -7);

  // Nodes within r cover the disc: about pi r^2 / step^2 of them.
  int count = 0;
  g.for_each_within(lat, lon, 1000.0, [&](uint64_t, double la, double lo) {
    CHECK(haversine_m(lat, lon, la, lo) <= 1000.0 * 1.01);
    ++count;
  });
  CHECK(count == Approx(3.14159 * 100).epsilon(0.1));
}
