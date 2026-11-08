// Doppler physics.
//
// Observer at rest, source moving through still air. The frequency heard at reception time t_r is
//     f(t_r) = f0 * c / (c - v_r(t_e)),      t_r = t_e + |p(t_e) - o| / c
// where v_r is the component of source velocity *toward* the observer, evaluated at the EMISSION
// time t_e (the "retarded" time) - not at t_r. Using v_r(t_r) produces a different, wrong curve.
#pragma once
#include <string>
#include <vector>

#include "adsb.h"
#include "dsp.h"
#include "geo.h"

namespace pl {

inline double speed_of_sound(double temp_c) { return 331.3 * std::sqrt(1.0 + temp_c / 273.15); }

// Straight-line constant-velocity flyover. v: speed, d: slant CPA distance, tc: audio time at which
// the CPA sound arrives (the "delayed CPA"). Returns ln(f/f0).
double straight_lng(double t_r, double v, double d, double tc, double c);

// Frequency ridge of one Doppler curve inside a user box.
struct Ridge {
  std::vector<double> t, f, w;  // audio time, Hz, weight (whitened brightness)
};
Ridge extract_ridge(const Spectrogram& whitened, double t0, double t1, double fmin, double fmax,
                    double max_slope_hz_s = 25.0, float min_level = 2.5f);

struct DopplerFit {
  double v = 0, d = 0, tc = 0, f0 = 0, cost = 0;
  double vmin = 0, vmax = 0, dmin = 0, dmax = 0, tcmin = 0, tcmax = 0, f0min = 0, f0max = 0;
  int npts = 0;
};
// side: 0 both, -1 only before CPA (approach), +1 only after CPA (recede).
// rel_tol: parameter sets whose robust cost <= best*(1+rel_tol)+abs_tol define the ranges.
DopplerFit fit_straight(const Ridge& r, double c, int side = 0, double rel_tol = 0.3, double abs_tol = 0.001);

// ln(f/f0) heard at observer o (ECEF) for absolute reception times trx[0..n). NaN when the emission
// time falls outside the track's coverage.
void model_lng(const Track& tr, const Vec3& o, const double* trx, int n, double c, float* out);

// Closest point of approach of a track to o, searched in [t_lo, t_hi]. Returns false unless a true
// interior local minimum exists. t_cpa: emission time, d: slant distance.
bool track_cpa(const Track& tr, const Vec3& o, double t_lo, double t_hi, double& t_cpa, double& d, double& speed);

}  // namespace pl
