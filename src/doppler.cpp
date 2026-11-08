#include "doppler.h"

#include <cmath>
#include <mutex>

#include "gpu.h"
#include "util.h"

namespace pl {

double straight_lng(double t_r, double v, double d, double tc, double c) {
  double t0 = tc - d / c;  // emission time of the CPA
  double T = t_r - t0;
  double c2 = c * c, a = c2 - v * v;
  double disc = c2 * v * v * T * T + a * d * d;
  double tau =
      (c2 * T - std::sqrt(disc)) / a;  // emission time relative to CPA emission
  double x = v * tau;
  double r = std::sqrt(x * x + d * d);
  double vr = -v * x / r;  // + when approaching
  return -std::log1p(-vr / c);
}

Ridge extract_ridge(const Spectrogram& s, double t0, double t1, double fmin,
                    double fmax, double max_slope, float min_level) {
  Ridge R;
  int f0 = std::max(0, (int)std::ceil((t0 - s.t0) / s.hop));
  int f1 = std::min(s.nframes - 1, (int)std::floor((t1 - s.t0) / s.hop));
  int b0 = std::max(0, (int)std::ceil(fmin / s.df) - s.bin0);
  int b1 = std::min(s.nbins - 1, (int)std::floor(fmax / s.df) - s.bin0);
  if (f1 <= f0 || b1 <= b0) return R;
  int nb = b1 - b0 + 1, nf = f1 - f0 + 1;
  int J = std::max(1, (int)std::ceil(max_slope * s.hop / s.df));
  const float lam = 0.15f;
  std::vector<float> acc(nb), nacc(nb);
  std::vector<int> back((size_t)nf * nb);
  for (int b = 0; b < nb; ++b)
    acc[b] = s.at(f0, b0 + b);
  for (int f = 1; f < nf; ++f) {
    for (int b = 0; b < nb; ++b) {
      float best = -1e30f;
      int arg = b;
      for (int k = std::max(0, b - J); k <= std::min(nb - 1, b + J); ++k) {
        float v = acc[k] - lam * std::abs(k - b);
        if (v > best) best = v, arg = k;
      }
      nacc[b] = best + s.at(f0 + f, b0 + b);
      back[(size_t)f * nb + b] = arg;
    }
    std::swap(acc, nacc);
  }
  int b = (int)(std::max_element(acc.begin(), acc.end()) - acc.begin());
  std::vector<int> path(nf);
  for (int f = nf - 1; f >= 0; --f) {
    path[f] = b;
    if (f) b = back[(size_t)f * nb + b];
  }
  for (int f = 0; f < nf; ++f) {
    int bb = b0 + path[f];
    float lv = s.at(f0 + f, bb);
    if (lv < min_level) continue;
    // parabolic interpolation of the peak in dB-z units
    double off = 0;
    if (bb > 0 && bb + 1 < s.nbins) {
      double ym = s.at(f0 + f, bb - 1), y0 = lv, yp = s.at(f0 + f, bb + 1);
      double den = ym - 2 * y0 + yp;
      if (den < 0) off = std::clamp(0.5 * (ym - yp) / den, -0.5, 0.5);
    }
    R.t.push_back(s.time(f0 + f));
    R.f.push_back(s.freq(bb) + off * s.df);
    R.w.push_back(lv);
  }
  return R;
}

namespace {
struct FitPts {
  std::vector<double> t, lnf, w;
};
// Robust cost of one parameter set; returns cost and the best ln f0.
inline double fit_cost(const FitPts& P, double v, double d, double tc, double c,
                       int side, double& lnf0) {
  const double cap = 0.02;
  size_t n = P.t.size();
  double sw = 0, s = 0;
  double res_buf[512];
  double* res = res_buf;
  std::vector<double> big;
  if (n > 512) big.resize(n), res = big.data();
  for (size_t i = 0; i < n; ++i) {
    if (side && (P.t[i] - tc) * side < 0) {
      res[i] = NAN;
      continue;
    }
    res[i] = P.lnf[i] - straight_lng(P.t[i], v, d, tc, c);
    s += P.w[i] * res[i];
    sw += P.w[i];
  }
  if (sw <= 0) return 1e9;
  double m = s / sw;
  for (int it = 0; it < 2; ++it) {
    double s2 = 0, w2 = 0;
    for (size_t i = 0; i < n; ++i)
      if (!std::isnan(res[i]) && std::fabs(res[i] - m) < cap)
        s2 += P.w[i] * res[i], w2 += P.w[i];
    if (w2 > 0) m = s2 / w2;
  }
  double cost = 0;
  for (size_t i = 0; i < n; ++i)
    if (!std::isnan(res[i]))
      cost += P.w[i] * std::min(std::fabs(res[i] - m), cap);
  lnf0 = m;
  return cost / sw;
}
}  // namespace

DopplerFit fit_straight(const Ridge& r, double c, int side, double rel_tol,
                        double abs_tol) {
  DopplerFit F;
  if (r.t.size() < 8)
    die("Doppler fit: too few ridge points (%zu) - widen the box or lower "
        "--min-level",
        r.t.size());
  FitPts P;
  size_t step = std::max<size_t>(1, r.t.size() / 300);
  for (size_t i = 0; i < r.t.size(); i += step) {
    P.t.push_back(r.t[i]);
    P.lnf.push_back(std::log(r.f[i]));
    P.w.push_back(r.w[i]);
  }
  F.npts = (int)P.t.size();
  double ta = P.t.front(), tb = P.t.back();

  struct Best {
    double cost = 1e30, v, d, tc, lnf0;
  };
  std::vector<float> gcost, glf;
  auto search = [&](std::vector<double> vs, std::vector<double> ds,
                    std::vector<double> tcs) {
    if (gpu::fit_costs(P.t, P.lnf, P.w, vs, ds, tcs, c, side, gcost, glf)) {
      Best b;
      size_t i = 0;
      for (double v : vs)
        for (double d : ds)
          for (double tc : tcs) {
            if (gcost[i] < b.cost) b = {gcost[i], v, d, tc, glf[i]};
            ++i;
          }
      return b;
    }
    int nt = num_threads();
    std::vector<Best> bests(nt);
    parallel_for(vs.size() * ds.size(), [&](size_t k, int tid) {
      double v = vs[k / ds.size()], d = ds[k % ds.size()];
      for (double tc : tcs) {
        double lf = 0;
        double cst = fit_cost(P, v, d, tc, c, side, lf);
        if (cst < bests[tid].cost) bests[tid] = {cst, v, d, tc, lf};
      }
    });
    Best b;
    for (auto& x : bests)
      if (x.cost < b.cost) b = x;
    return b;
  };
  auto lin = [](double a, double b, double st) {
    std::vector<double> o;
    for (double x = a; x <= b + 1e-9; x += st)
      o.push_back(x);
    return o;
  };
  auto geo = [](double a, double b, double fac) {
    std::vector<double> o;
    for (double x = a; x <= b * 1.0001; x *= fac)
      o.push_back(x);
    return o;
  };
  Best b = search(lin(15, 320, 2.5), geo(25, 8000, 1.06), lin(ta, tb, 0.5));
  b = search(lin(std::max(5.0, b.v - 5), b.v + 5, 0.25),
             geo(b.d * 0.85, b.d * 1.15, 1.005), lin(b.tc - 1, b.tc + 1, 0.05));
  F.v = b.v;
  F.d = b.d;
  F.tc = b.tc;
  F.f0 = std::exp(b.lnf0);
  F.cost = b.cost;

  // Uncertainty region: every parameter set whose cost is within tolerance of
  // the optimum.
  double thr = b.cost * (1 + rel_tol) + abs_tol;
  auto vs = lin(std::max(5.0, b.v * 0.5), b.v * 1.6, std::max(0.5, b.v * 0.01));
  auto ds = geo(std::max(10.0, b.d / 3), b.d * 3, 1.03);
  auto tcs = lin(b.tc - 6, b.tc + 6, 0.1);
  F.vmin = F.dmin = F.tcmin = F.f0min = 1e30;
  F.vmax = F.dmax = F.tcmax = F.f0max = -1e30;
  std::mutex mu;
  if (gpu::fit_costs(P.t, P.lnf, P.w, vs, ds, tcs, c, side, gcost, glf)) {
    size_t i = 0;
    for (double v : vs)
      for (double d : ds)
        for (double tc : tcs) {
          if (gcost[i] <= thr) {
            F.vmin = std::min(F.vmin, v), F.vmax = std::max(F.vmax, v);
            F.dmin = std::min(F.dmin, d), F.dmax = std::max(F.dmax, d);
            F.tcmin = std::min(F.tcmin, tc), F.tcmax = std::max(F.tcmax, tc);
            F.f0min = std::min(F.f0min, std::exp((double)glf[i])),
            F.f0max = std::max(F.f0max, std::exp((double)glf[i]));
          }
          ++i;
        }
    return F;
  }
  parallel_for(vs.size() * ds.size(), [&](size_t k, int) {
    double v = vs[k / ds.size()], d = ds[k % ds.size()];
    double lv[4] = {1e30, -1e30, 1e30, -1e30};  // tcmin tcmax f0min f0max
    bool any = false;
    for (double tc : tcs) {
      double lf = 0;
      if (fit_cost(P, v, d, tc, c, side, lf) <= thr) {
        any = true;
        lv[0] = std::min(lv[0], tc);
        lv[1] = std::max(lv[1], tc);
        lv[2] = std::min(lv[2], std::exp(lf));
        lv[3] = std::max(lv[3], std::exp(lf));
      }
    }
    if (!any) return;
    std::lock_guard<std::mutex> g(mu);
    F.vmin = std::min(F.vmin, v);
    F.vmax = std::max(F.vmax, v);
    F.dmin = std::min(F.dmin, d);
    F.dmax = std::max(F.dmax, d);
    F.tcmin = std::min(F.tcmin, lv[0]);
    F.tcmax = std::max(F.tcmax, lv[1]);
    F.f0min = std::min(F.f0min, lv[2]);
    F.f0max = std::max(F.f0max, lv[3]);
  });
  return F;
}

// Emission-sampled Doppler model. Instead of solving t_r = t_e + |p(t_e)-o|/c
// for t_e at every frame (a fixed-point iteration costing several trajectory
// evaluations per frame), sample the trajectory ONCE per step on a uniform
// EMISSION-time grid, compute when each sample's sound arrives (monotonic in
// t_e for subsonic sources) and interpolate ln(f/f0) onto the reception frames.
// Cost: ~1 trajectory evaluation per frame instead of ~4-10; error < 1e-4 in ln
// f at 0.1 s steps.
void model_lng(const Track& tr, const Vec3& o, const double* trx, int n,
               double c, float* out) {
  if (n <= 0) return;
  int hint = -1;
  Vec3 p, v;
  auto solve = [&](double t) {  // emission time for one reception time
    tr.eval(t, p, v, &hint);
    double te = t - (p - o).norm() / c;
    for (int k = 0; k < 8; ++k) {
      tr.eval(te, p, v, &hint);
      te = t - (p - o).norm() / c;
    }
    return te;
  };
  const double h = n > 1 ? (trx[n - 1] - trx[0]) / (n - 1) : 0.1;
  const double te_a = solve(trx[0]) - h, te_b = solve(trx[n - 1]) + h;
  const int ne = (int)std::ceil((te_b - te_a) / h) + 2;
  thread_local std::vector<double> tre;
  thread_local std::vector<float> lge;
  tre.resize(ne);
  lge.resize(ne);
  hint = -1;
  for (int j = 0; j < ne; ++j) {
    double te = te_a + j * h;
    tr.eval(te, p, v, &hint);
    Vec3 rv = o - p;
    double r = rv.norm();
    tre[j] = te + r / c;
    lge[j] = (te < tr.t0 || te > tr.t1)
                 ? NAN
                 : (float)(-std::log1p(-v.dot(rv) / r / c));
  }
  int j = 0;
  for (int i = 0; i < n; ++i) {
    double t = trx[i];
    while (j + 2 < ne && tre[j + 1] <= t)
      ++j;
    double w = (t - tre[j]) / (tre[j + 1] - tre[j]);
    if (w < 0 || w > 1 || std::isnan(lge[j]) || std::isnan(lge[j + 1])) {
      out[i] = NAN;
      continue;
    }
    out[i] = (float)(lge[j] + w * (lge[j + 1] - lge[j]));
  }
}

bool track_cpa(const Track& tr, const Vec3& o, double t_lo, double t_hi,
               double& t_cpa, double& d, double& speed) {
  double a = std::max(t_lo, tr.t0), b = std::min(t_hi, tr.t1);
  if (b - a < 2) return false;
  int hint = -1;
  Vec3 p, v;
  double best = 1e30, bt = a;
  int n = (int)std::ceil(b - a) + 1;
  int bi = 0;
  for (int i = 0; i < n; ++i) {
    double t = std::min(b, a + i);
    tr.eval(t, p, v, &hint);
    double dd = (p - o).norm2();
    if (dd < best) best = dd, bt = t, bi = i;
  }
  if (bi == 0 || bi == n - 1) return false;
  double lo = bt - 1, hi = bt + 1;
  // The Hermite velocity is the exact derivative of the position, so the CPA is
  // the root of f(t) = v(t).(p(t) - o) (half the derivative of the squared
  // range). Illinois regula falsi converges in ~5 evaluations to the same
  // tolerance golden section needs 32 for; golden section stays as fallback.
  {
    auto F = [&](double t) {
      tr.eval(t, p, v, &hint);
      return v.dot(p - o);
    };
    double flo = F(lo), fhi = F(hi);
    if (flo < 0 && fhi > 0) {
      double t = bt, tp = 1e300;
      int side = 0;
      for (int it = 0; it < 40 && hi - lo > 1e-6 && std::fabs(t - tp) > 1e-9;
           ++it) {
        tp = t;
        t = (lo * fhi - hi * flo) / (fhi - flo);
        double ft = F(t);
        if (ft == 0) break;
        if (ft < 0) {
          lo = t, flo = ft;
          if (side == -1) fhi *= 0.5;
          side = -1;
        } else {
          hi = t, fhi = ft;
          if (side == 1) flo *= 0.5;
          side = 1;
        }
      }
      t_cpa = t;
      tr.eval(t_cpa, p, v, &hint);
      d = (p - o).norm();
      speed = v.norm();
      return true;
    }
  }
  // golden section on [bt-1, bt+1]
  const double g = 0.618033988749895;
  double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
  auto D = [&](double t) {
    tr.eval(t, p, v, &hint);
    return (p - o).norm2();
  };
  double f1 = D(x1), f2 = D(x2);
  for (int it = 0; it < 30; ++it) {
    if (f1 < f2)
      hi = x2, x2 = x1, f2 = f1, x1 = hi - g * (hi - lo), f1 = D(x1);
    else
      lo = x1, x1 = x2, f1 = f2, x2 = lo + g * (hi - lo), f2 = D(x2);
  }
  t_cpa = 0.5 * (lo + hi);
  tr.eval(t_cpa, p, v, &hint);
  d = (p - o).norm();
  speed = v.norm();
  return true;
}

}  // namespace pl
