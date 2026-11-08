#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "doppler.h"

using Catch::Approx;
using namespace pl;

TEST_CASE("speed of sound follows temperature", "[doppler]") {
  CHECK(speed_of_sound(0) == Approx(331.3));
  CHECK(speed_of_sound(20) == Approx(343.2).margin(0.1));
}

TEST_CASE("straight-line Doppler curve", "[doppler]") {
  const double c = 343, v = 60, d = 900, tc = 300;
  // Approaching is higher, receding lower, and the far limits are the
  // classical c/(c -+ v) shifts.
  CHECK(straight_lng(tc - 200, v, d, tc, c) > 0);
  CHECK(straight_lng(tc + 200, v, d, tc, c) < 0);
  CHECK(straight_lng(tc - 1e6, v, d, tc, c) ==
        Approx(std::log(c / (c - v))).margin(1e-6));
  CHECK(straight_lng(tc + 1e6, v, d, tc, c) ==
        Approx(std::log(c / (c + v))).margin(1e-6));
  // tc is the arrival of the sound emitted at closest approach, where the
  // radial velocity is zero.
  CHECK(straight_lng(tc, v, d, tc, c) == Approx(0.0).margin(1e-9));
  // A closer pass sweeps faster through the centre.
  const double near = straight_lng(tc + 5, v, d / 3, tc, c);
  const double far = straight_lng(tc + 5, v, d, tc, c);
  CHECK(near < far);
}

TEST_CASE("fit_straight recovers a synthetic flyover", "[doppler]") {
  const double c = 343, v = 55, d = 850, tc = 120, f0 = 110;
  Ridge r;
  for (double t = 40; t <= 200; t += 0.25) {
    r.t.push_back(t);
    r.f.push_back(f0 * std::exp(straight_lng(t, v, d, tc, c)));
    r.w.push_back(5.0);
  }
  const DopplerFit fit = fit_straight(r, c);
  CHECK(fit.v == Approx(v).epsilon(0.03));
  CHECK(fit.d == Approx(d).epsilon(0.05));
  CHECK(fit.tc == Approx(tc).margin(0.5));
  CHECK(fit.f0 == Approx(f0).epsilon(0.005));
  CHECK(fit.vmin <= fit.v);
  CHECK(fit.v <= fit.vmax);
  CHECK(fit.dmin <= fit.d);
  CHECK(fit.d <= fit.dmax);
}
