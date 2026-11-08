// CUDA implementation of the hot paths. Each kernel mirrors the CPU algorithm
// (see the referenced functions) so results agree to rounding; the CPU code
// stays the reference implementation.
#include <climits>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "gpu.h"
#include "util.h"

namespace pl::gpu {
namespace {

#define CK(x)                                             \
  do {                                                    \
    cudaError_t e_ = (x);                                 \
    if (e_ != cudaSuccess) {                              \
      log("gpu: %s failed: %s (falling back to CPU)", #x, \
          cudaGetErrorString(e_));                        \
      g_enabled = false;                                  \
      return fail_value;                                  \
    }                                                     \
  } while (0)

bool g_enabled = false;
std::string g_info = "no CUDA device";

constexpr double kPiD = 3.14159265358979323846;
constexpr double kDegD = kPiD / 180.0;
constexpr double kA = 6378137.0, kF = 1.0 / 298.257223563,
                 kE2 = kF * (2.0 - kF);

// --------------------------------------------------------- device data ---
struct DTracks {
  const double* t = nullptr;
  const double* p = nullptr;  // xyz interleaved
  const double* v = nullptr;
  const uint64_t* off = nullptr;  // ntracks+1
  const double* t0 = nullptr;
  const double* t1 = nullptr;
  uint32_t n = 0;
};
struct DDem {
  const int* idx = nullptr;  // (la - la0) * nlo + (lo - lo0) -> tile or -1
  const int16_t* data = nullptr;
  const uint64_t* toff = nullptr;
  const int* tn = nullptr;
  int la0 = 0, lo0 = 0, nla = 0, nlo = 0;
  double def = 0;
  bool any = false;
};
struct DEvent {
  const float* himg;  // harmonic-summed sharp table [(j*nb + b)*nrH + rr]  (j =
                      // frame - fb0)
  const float* hblur3;  // harmonic-summed blur, phase split + zero padded:
                        // [((j*nb + b)*3 + p)*L3 + Z + q]
  int nb, E, nrH, nq, L3, Z;
  float thr[8];
  int nf, nr, fb0, fb1, rlo, rhi, nm;
  float mrow[8];
  double t0, hop, lnf0, dln, tc;
  double wpre, wpost;  // trajectory window (s) around the event's CPA guess
                       // needed by CPA search + model
};

DTracks g_T;
DDem g_D;
DEvent g_E[kMaxEvents];
int g_K = 0, g_S = 0;
double g_hop = 0.1, g_c = 343, g_obs = 2;
std::vector<void*> g_track_bufs, g_dem_bufs, g_event_bufs;
size_t g_max_rows = 0;  // max(nfr + 2S) over events (shared memory sizing)
int g_ne_cap = 0;       // emission-sample scratch capacity

template <class T>
void free_all(std::vector<T*>& v) {
  for (auto* p : v)
    cudaFree(p);
  v.clear();
}
void free_all(std::vector<void*>& v) {
  for (auto* p : v)
    cudaFree(p);
  v.clear();
}

// ------------------------------------------------------ device helpers ---
struct D3 {
  double x, y, z;
};
__device__ __forceinline__ D3 operator+(D3 a, D3 b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
__device__ __forceinline__ D3 operator-(D3 a, D3 b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
__device__ __forceinline__ D3 operator*(D3 a, double s) {
  return {a.x * s, a.y * s, a.z * s};
}
__device__ __forceinline__ double dot(D3 a, D3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
__device__ __forceinline__ double norm(D3 a) {
  return sqrt(dot(a, a));
}

__device__ D3 ld3(const double* a, uint64_t i) {
  return {a[3 * i], a[3 * i + 1], a[3 * i + 2]};
}

// Track::eval
__device__ void eval(const DTracks& T, uint32_t tr, double t, D3& p, D3& v,
                     int& hint) {
  uint64_t o0 = T.off[tr];
  int n = (int)(T.off[tr + 1] - o0);
  const double* tt = T.t + o0;
  if (t <= tt[0]) {
    D3 p0 = ld3(T.p, o0), v0 = ld3(T.v, o0);
    p = p0 + v0 * (t - tt[0]);
    v = v0;
    return;
  }
  if (t >= tt[n - 1]) {
    D3 p0 = ld3(T.p, o0 + n - 1), v0 = ld3(T.v, o0 + n - 1);
    p = p0 + v0 * (t - tt[n - 1]);
    v = v0;
    return;
  }
  int i;
  if (hint >= 0 && hint < n - 1) {
    i = hint;
    while (i < n - 2 && tt[i + 1] <= t)
      ++i;
    while (i > 0 && tt[i] > t)
      --i;
  } else {
    int lo = 0, hi = n;  // upper_bound
    while (lo < hi) {
      int mid = (lo + hi) >> 1;
      if (t < tt[mid])
        hi = mid;
      else
        lo = mid + 1;
    }
    i = lo - 1;
  }
  hint = i;
  D3 pa = ld3(T.p, o0 + i), va = ld3(T.v, o0 + i), pb = ld3(T.p, o0 + i + 1),
     vb = ld3(T.v, o0 + i + 1);
  double h = tt[i + 1] - tt[i], s = (t - tt[i]) / h, s2 = s * s, s3 = s2 * s;
  double h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s,
         h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
  double d00 = 6 * s2 - 6 * s, d10 = 3 * s2 - 4 * s + 1, d01 = -6 * s2 + 6 * s,
         d11 = 3 * s2 - 2 * s;
  p = pa * h00 + va * (h10 * h) + pb * h01 + vb * (h11 * h);
  v = (pa * d00 + va * (d10 * h) + pb * d01 + vb * (d11 * h)) * (1.0 / h);
}

__device__ D3 geo2ecef(double lat, double lon, double h) {
  double la = lat * kDegD, lo = lon * kDegD;
  double sl = sin(la), cl = cos(la);
  double N = kA / sqrt(1.0 - kE2 * sl * sl);
  return {(N + h) * cl * cos(lo), (N + h) * cl * sin(lo),
          (N * (1.0 - kE2) + h) * sl};
}

__device__ void ecef2geo(D3 p, double& lat_deg, double& lon_deg, double& h) {
  double b = kA * (1.0 - kF);
  double ep2 = (kA * kA - b * b) / (b * b);
  double r = sqrt(p.x * p.x + p.y * p.y);
  double lon = atan2(p.y, p.x);
  double th = atan2(p.z * kA, r * b);
  double st = sin(th), ct = cos(th);
  double lat = atan2(p.z + ep2 * b * st * st * st, r - kE2 * kA * ct * ct * ct);
  for (int it = 0; it < 2; ++it) {
    double sl = sin(lat);
    double N = kA / sqrt(1.0 - kE2 * sl * sl);
    double hh = r / cos(lat) - N;
    lat = atan2(p.z, r * (1.0 - kE2 * N / (N + hh)));
  }
  double sl = sin(lat);
  double N = kA / sqrt(1.0 - kE2 * sl * sl);
  h = (fabs(cos(lat)) > 1e-9) ? r / cos(lat) - N : fabs(p.z) - b;
  lat_deg = lat / kDegD;
  lon_deg = lon / kDegD;
}

__device__ double mlat_of(double lat) {
  double s = sin(lat * kDegD);
  double w = 1.0 - kE2 * s * s;
  return kDegD * kA * (1.0 - kE2) / (w * sqrt(w));
}
__device__ double mlon_of(double lat) {
  double s = sin(lat * kDegD);
  return kDegD * kA * cos(lat * kDegD) / sqrt(1.0 - kE2 * s * s);
}
__device__ double row_dlon(double step, double dlat, int i) {
  double c = cos(i * dlat * kDegD);
  if (c < 1e-6) c = 1e-6;
  return step / (111319.49 * c);
}
__device__ __forceinline__ uint64_t pack(int32_t i, int32_t j) {
  return ((uint64_t)(uint32_t)i << 32) | (uint32_t)j;
}

// Dem::elev
__device__ double elev(const DDem& D, double lat, double lon) {
  if (!D.any) return D.def;
  int la = (int)floor(lat), lo = (int)floor(lon);
  int a = la - D.la0, b = lo - D.lo0;
  if (a < 0 || b < 0 || a >= D.nla || b >= D.nlo) return D.def;
  int id = D.idx[a * D.nlo + b];
  if (id < 0) return D.def;
  int n = D.tn[id];
  const int16_t* hgt = D.data + D.toff[id];
  double x = (lon - lo) * (n - 1), y = (la + 1 - lat) * (n - 1);
  int ix = min(max((int)x, 0), n - 2), iy = min(max((int)y, 0), n - 2);
  double fx = x - ix, fy = y - iy;
  auto H = [&](int xx, int yy) {
    int16_t v = hgt[(size_t)yy * n + xx];
    return v == -32768 ? D.def : (double)v;
  };
  return (H(ix, iy) * (1 - fx) + H(ix + 1, iy) * fx) * (1 - fy) +
         (H(ix, iy + 1) * (1 - fx) + H(ix + 1, iy + 1) * fx) * fy;
}

// track_cpa (serial; used per candidate in stage A)
__device__ bool cpa_serial(const DTracks& T, uint32_t tr, D3 o, double t_lo,
                           double t_hi, double& t_cpa, double& d,
                           double& speed) {
  double a = fmax(t_lo, T.t0[tr]), b = fmin(t_hi, T.t1[tr]);
  if (b - a < 2) return false;
  int hint = -1;
  D3 p, v;
  double best = 1e30, bt = a;
  int n = (int)ceil(b - a) + 1, bi = 0;
  for (int i = 0; i < n; ++i) {
    double t = fmin(b, a + i);
    eval(T, tr, t, p, v, hint);
    D3 r = p - o;
    double dd = dot(r, r);
    if (dd < best) best = dd, bt = t, bi = i;
  }
  if (bi == 0 || bi == n - 1) return false;
  double lo = bt - 1, hi = bt + 1;
  {  // Illinois regula falsi on f = v.(p - o), as track_cpa on the CPU
    auto F = [&](double t) {
      eval(T, tr, t, p, v, hint);
      return dot(v, p - o);
    };
    double flo = F(lo), fhi = F(hi);
    if (flo < 0 && fhi > 0) {
      double t = bt, tp = 1e300;
      int side = 0;
      for (int it = 0; it < 40 && hi - lo > 1e-6 && fabs(t - tp) > 1e-9; ++it) {
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
      eval(T, tr, t_cpa, p, v, hint);
      d = norm(p - o);
      speed = norm(v);
      return true;
    }
  }
  const double g = 0.618033988749895;
  double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
  auto Dd = [&](double t) {
    eval(T, tr, t, p, v, hint);
    D3 r = p - o;
    return dot(r, r);
  };
  double f1 = Dd(x1), f2 = Dd(x2);
  for (int it = 0; it < 30; ++it) {
    if (f1 < f2)
      hi = x2, x2 = x1, f2 = f1, x1 = hi - g * (hi - lo), f1 = Dd(x1);
    else
      lo = x1, x1 = x2, f1 = f2, x2 = lo + g * (hi - lo), f2 = Dd(x2);
  }
  t_cpa = 0.5 * (lo + hi);
  eval(T, tr, t_cpa, p, v, hint);
  d = norm(p - o);
  speed = norm(v);
  return true;
}

// straight_lng
__device__ double straight_lng(double t_r, double v, double d, double tc,
                               double c) {
  double t0 = tc - d / c, T = t_r - t0;
  double c2 = c * c, a = c2 - v * v;
  double disc = c2 * v * v * T * T + a * d * d;
  double tau = (c2 * T - sqrt(disc)) / a;
  double x = v * tau, r = sqrt(x * x + d * d);
  double vr = -v * x / r;
  return -log1p(-vr / c);
}

// ======================================================= stage A kernels ===
struct Samp {
  double lat, lon, h, mlat, mlon;
  uint32_t range;
  int32_t k;
  int valid;
};

__global__ void k_samples(DTracks T, DDem D, const SampleRange* ranges,
                          const uint64_t* prefix, int nr, uint64_t N,
                          HitParams hp, double Rh, double mla, double mlo,
                          Samp* out) {
  uint64_t s = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  if (s >= N) return;
  int lo = 0, hi = nr;  // upper_bound(prefix, s) - 1
  while (lo < hi) {
    int mid = (lo + hi) >> 1;
    if (s < prefix[mid])
      hi = mid;
    else
      lo = mid + 1;
  }
  int r = lo - 1;
  const SampleRange R = ranges[r];
  int k = R.ka + (int)(s - prefix[r]);
  Samp o{};
  o.range = r;
  o.k = k;
  o.valid = 0;
  double t = R.a + k * hp.dt;
  int hint = -1;
  D3 p, v;
  eval(T, R.track, t, p, v, hint);
  double sp = norm(v);
  if (!(sp < hp.vmin * 0.97 || sp > hp.vmax * 1.03)) {
    double lat, lon, h;
    ecef2geo(p, lat, lon, h);
    if (!(lat < hp.la0 - mla || lat > hp.la1 + mla || lon < hp.lo0 - mlo ||
          lon > hp.lo1 + mlo) &&
        !(h - elev(D, lat, lon) - 600 > Rh)) {
      o.lat = lat, o.lon = lon, o.h = h, o.mlat = mlat_of(lat),
      o.mlon = mlon_of(lat);
      o.valid = 1;
    }
  }
  out[s] = o;
}

struct Cand {
  uint64_t key;
  uint64_t s;
};

__device__ __forceinline__ float dist_or_inf(const Samp& q, double la,
                                             double lo, double ge, double obs,
                                             double Rh) {
  double dz = q.h - (ge + obs);
  double dx = (lo - q.lon) * q.mlon, dy = (la - q.lat) * q.mlat;
  double d2 = dx * dx + dy * dy + dz * dz;
  return d2 <= Rh * Rh ? (float)sqrt(d2) : INFINITY;
}

// One thread per (sample, grid row of its disk). Local-minimum test against the
// neighbouring samples replaces the CPU's sort-and-scan (identical semantics:
// absent neighbour == +inf).
__global__ void k_raster(DDem D, const Samp* S, uint64_t N, int maxRows,
                         HitParams hp, double Rh, double dlat, Cand* cand,
                         unsigned long long* ncand, unsigned long long cap) {
  uint64_t idx = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  uint64_t s = idx / maxRows;
  int rr = (int)(idx % maxRows);
  if (s >= N) return;
  const Samp q = S[s];
  if (!q.valid) return;
  double mlat = mlat_of(q.lat);
  int32_t i0 = (int32_t)ceil((q.lat - Rh / mlat) / dlat);
  int32_t i1 = (int32_t)floor((q.lat + Rh / mlat) / dlat);
  int32_t i = i0 + rr;
  if (i > i1) return;
  double la = i * dlat;
  double dy = (la - q.lat) * mlat;
  double rem = Rh * Rh - dy * dy;
  if (rem < 0) return;
  double dxm = sqrt(rem), mlonr = mlon_of(la), dl = row_dlon(hp.step, dlat, i);
  int32_t j0 = (int32_t)ceil((q.lon - dxm / mlonr) / dl),
          j1 = (int32_t)floor((q.lon + dxm / mlonr) / dl);
  bool hasPrev = s > 0 && S[s - 1].valid && S[s - 1].range == q.range &&
                 S[s - 1].k == q.k - 1;
  bool hasNext = s + 1 < N && S[s + 1].valid && S[s + 1].range == q.range &&
                 S[s + 1].k == q.k + 1;
  Samp qp, qn;
  if (hasPrev) qp = S[s - 1];
  if (hasNext) qn = S[s + 1];
  for (int32_t j = j0; j <= j1; ++j) {
    double lo = j * dl;
    if (!(la >= hp.la0 && la < hp.la1 && lo >= hp.lo0 && lo < hp.lo1)) continue;
    double ge = elev(D, la, lo);
    float d = dist_or_inf(q, la, lo, ge, hp.obs_h, Rh);
    if (isinf(d)) continue;
    float dp = hasPrev ? dist_or_inf(qp, la, lo, ge, hp.obs_h, Rh) : INFINITY;
    if (!(dp >= d)) continue;
    float dn = hasNext ? dist_or_inf(qn, la, lo, ge, hp.obs_h, Rh) : INFINITY;
    if (!(dn > d)) continue;
    unsigned long long w = atomicAdd(ncand, 1ULL);
    if (w < cap) cand[w] = {pack(i, j), s};
  }
}

__global__ void k_refine(DTracks T, DDem D, const Samp* S,
                         const SampleRange* ranges, const Cand* cand,
                         uint64_t nc, HitParams hp, double dlat, GHit* hits,
                         unsigned long long* nh, unsigned long long cap) {
  uint64_t x = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
  if (x >= nc) return;
  Cand cd = cand[x];
  const Samp& q = S[cd.s];
  const SampleRange R = ranges[q.range];
  int32_t i = (int32_t)(uint32_t)(cd.key >> 32),
          j = (int32_t)(uint32_t)(cd.key & 0xffffffffu);
  double la = i * dlat, lo = j * row_dlon(hp.step, dlat, i);
  D3 o = geo2ecef(la, lo, elev(D, la, lo) + hp.obs_h);
  double tk = R.a + q.k * hp.dt, tcpa, d, sp;
  if (!cpa_serial(T, R.track, o, tk - 2 * hp.dt - 1, tk + 2 * hp.dt + 1, tcpa,
                  d, sp))
    return;
  if (d < hp.dmin || d > hp.dmax || sp < hp.vmin || sp > hp.vmax) return;
  double Tt = tcpa + d / hp.c;
  if (Tt < hp.Tlo || Tt > hp.Thi) return;
  unsigned long long w = atomicAdd(nh, 1ULL);
  if (w < cap) hits[w] = {cd.key, R.track, (float)d, (float)sp, Tt};
}

// ======================================================= stage D kernel ===
constexpr int kBlock = 256;
constexpr int kWarps = kBlock / 32;
constexpr int kCombo = 1024;
constexpr int kGS =
    16;  // shifts per coarse alignment pass (register accumulators)

__constant__ DEvent c_E[kMaxEvents];

#ifdef PL_PROF_PHASES
__device__ unsigned long long
    g_phase[12];  // 0 load, 1 cpa, 2 model, 3 align, 4 fallback, 5 jobs, 6
                  // coarse-eval, 7 coarse-serial, 8 refine-eval, 9
                  // refine-serial
#define PHASE_T0() long long _pt = clock64()
#define PHASE_ADD(i)                                                 \
  do {                                                               \
    if (threadIdx.x == 0)                                            \
      atomicAdd(&g_phase[i], (unsigned long long)(clock64() - _pt)); \
    _pt = clock64();                                                 \
  } while (0)
#define PHASE_CNT(i)                                    \
  do {                                                  \
    if (threadIdx.x == 0) atomicAdd(&g_phase[i], 1ULL); \
  } while (0)
#else
#define PHASE_T0()
#define PHASE_ADD(i)
#define PHASE_CNT(i)
#endif

struct HitLess {
  __host__ __device__ bool operator()(const GHit& a, const GHit& b) const {
    return a.key != b.key       ? a.key < b.key
           : a.track != b.track ? a.track < b.track
                                : a.T < b.T;
  }
};
struct HitDup {  // same node + flight within 2 s: one pass reported twice
                 // (range split / overlap)
  __host__ __device__ bool operator()(const GHit& a, const GHit& b) const {
    return a.key == b.key && a.track == b.track && fabs(b.T - a.T) < 2.0;
  }
};

struct Top {
  float sc;
  int s, r;
};
// Insert into a top-3 list ordered by (score desc, index asc): the same set and
// order a sequential "replace the third if strictly greater" scan in index
// order produces.
__device__ __forceinline__ void top3_push(Top* t, Top x) {
  auto better = [](const Top& a, const Top& b) {
    return a.sc > b.sc || (a.sc == b.sc && a.s < b.s);
  };
  if (!better(x, t[2])) return;
  if (better(x, t[1])) {
    t[2] = t[1];
    if (better(x, t[0]))
      t[1] = t[0], t[0] = x;
    else
      t[1] = x;
  } else
    t[2] = x;
}

__device__ __forceinline__ int harm_variant(const DEvent& e, float x, int& fl) {
  float f = floorf(x);
  fl = (int)f;
  float phi = x - f;
  int b = 0;
  while (b < e.nb - 1 && phi >= e.thr[b])
    ++b;
  return b;
}

// Coarse sweep inner loop: frames [ja, jb) (coarse-frame units) x kGS shifts,
// one harmonic-table read each. Loads are unconditional (no-model rows point
// into zero padding; surplus shifts are clamped duplicates whose sums are
// discarded) so the unrolled body issues its loads back to back.
__device__ __forceinline__ void coarse_frames(const DEvent& e, const int* bo,
                                              int ja, int jb, int S, int g0,
                                              int gn, int ss, int k,
                                              float* acc) {
  const size_t fstride = (size_t)e.nb * 3 * e.L3;
  for (int jj = ja; jj < jb; ++jj) {
    const int j = jj * 3;
    const float* P = e.hblur3 + (size_t)j * fstride + k;
    const int* B = bo + (j + 2 * S - g0 * ss);
#pragma unroll
    for (int i = 0; i < kGS; ++i)
      acc[i] += __ldg(P + B[-min(i, gn - 1) * ss]);
  }
}

// Refine inner loop: frames ja, ja+4, .. < jb x 15 (seed, shift) segments, lane
// r offset ri. rf/rv: per model row floor row + E (INT_MIN/4 = no model) and
// variant base b*nrH, precomputed once per job.
__device__ __forceinline__ void refine_frames(const DEvent& e, const int* rf,
                                              const int* rv, int ja, int jb,
                                              int S, int nseg, const int* seg_s,
                                              const int* seg_r, int ri,
                                              float* acc) {
  const int nrH = e.nrH;
  for (int j = ja; j < jb; j += 4) {
    const float* P = e.himg + (size_t)j * e.nb * nrH;
#pragma unroll
    for (int sg = 0; sg < 15; ++sg) {
      const int g = min(sg, nseg - 1);
      const int i = j + S - seg_s[g];
      const int row = rf[i] + seg_r[g] + ri;
      const int rc = min(max(row, 0), nrH - 1);
      const float v = __ldg(P + rv[i] + rc);
      acc[sg] += row == rc ? v : 0.f;
    }
  }
}

// Alignment sweep (see best_alignment in search.cpp). lrint(x + m + r) ==
// lrint(x + m) + r, so for a fixed shift s every f0 offset r shares one base
// row per (frame, harmonic): warp lanes take consecutive r and read consecutive
// floats of one image row (phase-split blur for the step-3 coarse grid) -> one
// 128-byte line per warp load instead of ~32 scattered sectors. Each warp
// accumulates a frame chunk in registers and folds it into the shared combo
// scores with one shared atomic per lane.
__device__ float block_align(const DEvent& e, const float* rows, int nfr, int S,
                             int& bs, int& br, float* scratch, int scratchN) {
  const int warp = threadIdx.x >> 5, lane = threadIdx.x & 31;
  __shared__ Top top[3];
  __shared__ float s_best;
  __shared__ int s_bs, s_br;
  __shared__ int seg_s[15], seg_r[15], seg_w[15], s_nseg;
  const int ss = S >= 2 ? 2 : 1, nm = e.nm, nq = e.nq;
  int ns = 0;
  for (int s = -S; s <= S; s += ss)
    ++ns;
  const int nrr = e.rhi >= e.rlo ? (e.rhi - e.rlo) / 3 + 1 : 0;
  const int nkb = (nrr + 31) >> 5;
  const int ncf = (nfr + 2) / 3;  // coarse frames
  const float norm_c = (float)(ncf * nm);
  // Per-job precompute, shared by all 32*nkb lanes and all shifts: for every
  // model row i and harmonic m the offset of r = rlo inside a frame of the
  // padded phase-split blur (-1 = no model). Turns the inner loop into
  // shared-load + global-load + add. Lives in the (now idle)
  // trajectory/emission scratch.
  const int nrow = nfr + 2 * S;
  int* bo = reinterpret_cast<int*>(scratch);
  float* sc = scratch + nrow;
  const int scCap = scratchN - nrow;
  for (int i = threadIdx.x; i < nrow; i += blockDim.x) {
    float x = rows[i];
    int v = 0;  // no model: offset 0 = left zero padding of (variant 0, phase
                // 0): reads 0 for every lane
    if (!isnan(x)) {
      int fl, b = harm_variant(e, x, fl);
      int row0 = fl + e.rlo + e.E;
      int p = ((row0 % 3) + 3) % 3, q0 = min(max((row0 - p) / 3, -e.Z), nq);
      v = (b * 3 + p) * e.L3 + e.Z + q0;
    }
    bo[i] = v;
  }
  __shared__ Top wtop[3 * kWarps];
  if (threadIdx.x == 0)
    for (int q = 0; q < 3; ++q)
      top[q] = {-1.f, -1, 0};
  __syncthreads();
  PHASE_T0();
  // Frame-outer / shift-inner: a warp owns 32 consecutive r and a frame chunk,
  // and for every frame visits all shifts of the group, so one image frame is
  // pulled from L2 once and then re-read from L1 by every shift (~31x less L2
  // traffic than combo-outer). Per-shift partial sums live in registers.
  const int gs = max(1, min(min(ns, kGS), scCap / max(nrr, 1)));
  for (int g0 = 0; g0 < ns; g0 += gs) {
    const int gn = min(gs, ns - g0), C = gn * nrr;
    for (int c = threadIdx.x; c < C; c += blockDim.x)
      sc[c] = 0.f;
    __syncthreads();
    const int nch = max(1, min(ncf / 4, (3 * kWarps + nkb - 1) / nkb));
    const int ntask = nkb * nch;
    for (int t = warp; t < ntask; t += kWarps) {
      int ch = t % nch, kb = t / nch;
      int k = kb * 32 + lane;
      int ja = (int)((long long)ncf * ch / nch),
          jb = (int)((long long)ncf * (ch + 1) / nch);
      float acc[kGS];
#pragma unroll
      for (int i = 0; i < kGS; ++i)
        acc[i] = 0.f;
      coarse_frames(e, bo, ja, jb, S, g0, gn, ss, k, acc);
      if (k < nrr)
#pragma unroll
        for (int i = 0; i < kGS; ++i)
          if (i < gn) atomicAdd(&sc[i * nrr + k], acc[i]);
    }
    __syncthreads();
    PHASE_ADD(6);
    {  // parallel top-3 over this group: per-lane lists -> warp shuffle merge
       // -> warp lists -> thread 0
      Top l[3] = {{-1.f, -1, 0},
                  {-1.f, -1, 0},
                  {-1.f, -1, 0}};  // .s = global combo index while ranking
      for (int c = threadIdx.x; c < C; c += blockDim.x)
        top3_push(l, {sc[c] / norm_c, g0 * nrr + c, 0});
      for (int o = 16; o; o >>= 1) {
        Top t[3];
        for (int q = 0; q < 3; ++q)
          t[q] = {__shfl_down_sync(0xffffffffu, l[q].sc, o),
                  __shfl_down_sync(0xffffffffu, l[q].s, o), 0};
        if (lane < o)
          for (int q = 0; q < 3; ++q)
            top3_push(l, t[q]);
      }
      if (lane == 0)
        for (int q = 0; q < 3; ++q)
          wtop[warp * 3 + q] = l[q];
      __syncthreads();
      if (threadIdx.x == 0)
        for (int q = 0; q < 3 * kWarps; ++q)
          top3_push(top, wtop[q]);
    }
    __syncthreads();
    PHASE_ADD(7);
  }
  // refine on the sharp image: segments (seed, s) x up to 7 consecutive r. Lane
  // = (frame phase 0..3, r 0..7).
  if (threadIdx.x == 0) {
    for (int q = 0; q < 3; ++q) {  // combo index -> (s, r)
      int c = top[q].s;
      top[q].s = c >= 0 ? -S + (c / nrr) * ss : 0;
      top[q].r = c >= 0 ? e.rlo + (c % nrr) * 3 : e.rlo;
    }
    int n = 0;
    for (int q = 0; q < 3; ++q) {
      if (top[q].sc < 0) continue;
      for (int s = max(-S, top[q].s - ss); s <= min(S, top[q].s + ss); ++s) {
        int ra = max(e.rlo, top[q].r - 3), rb = min(e.rhi, top[q].r + 3);
        seg_s[n] = s, seg_r[n] = ra, seg_w[n] = rb - ra + 1, ++n;
      }
    }
    s_nseg = n;
  }
  __syncthreads();
  const int nseg = s_nseg;
  for (int c = threadIdx.x; c < nseg * 8; c += blockDim.x)
    sc[c] = 0.f;
  // per-row refine offsets (reuse the coarse offset scratch: bo is dead now;
  // sc[0..120) sits after it)
  int* rf = bo;
  int* rv = reinterpret_cast<int*>(sc) + 128;
  for (int i = threadIdx.x; i < nrow; i += blockDim.x) {
    float x = rows[i];
    int fl = 0, b = 0;
    if (!isnan(x))
      b = harm_variant(e, x, fl), fl += e.E;
    else
      fl = INT_MIN / 4;
    rf[i] = fl, rv[i] = b * e.nrH;
  }
  __syncthreads();
  PHASE_ADD(9);
  {  // frame-outer / segment-inner (L1 reuse of each frame across the <=15
     // segments), lane = (frame 0..3, r 0..7)
    const int nch = max(1, min(nfr / 16, 3 * kWarps));
    const int fi = lane >> 3, ri = lane & 7;
    for (int ch = warp; ch < nch; ch += kWarps) {
      int ja = (int)((long long)nfr * ch / nch),
          jb = (int)((long long)nfr * (ch + 1) / nch);
      float acc[15];
#pragma unroll
      for (int i = 0; i < 15; ++i)
        acc[i] = 0.f;
      refine_frames(e, rf, rv, ja + fi, jb, S, nseg, seg_s, seg_r, ri, acc);
#pragma unroll
      for (int sg = 0; sg < 15; ++sg) {
        float a = acc[sg];
        a += __shfl_xor_sync(0xffffffffu, a, 8);
        a += __shfl_xor_sync(0xffffffffu, a, 16);
        if (fi == 0 && sg < nseg && ri < seg_w[sg])
          atomicAdd(&sc[sg * 8 + ri], a);
      }
    }
  }
  __syncthreads();
  PHASE_ADD(8);
  if (threadIdx.x == 0) {
    float best = -1;
    int b_s = 0, b_r = e.rlo;
    const float norm_f = (float)(nfr * nm);
    for (int sg = 0; sg < nseg; ++sg)
      for (int i = 0; i < seg_w[sg]; ++i) {
        float v = sc[sg * 8 + i] / norm_f;
        if (v > best) best = v, b_s = seg_s[sg], b_r = seg_r[sg] + i;
      }
    s_best = best, s_bs = b_s, s_br = b_r;
  }
  __syncthreads();
  PHASE_ADD(9);
  bs = s_bs;
  br = s_br;
  return s_best;
}

// model_lng for one frame (independent start; converges to the same emission
// time as the CPU loop)
__device__ float frame_lng(const DTracks& T, uint32_t tr, D3 o, double t,
                           double c) {
  int hint = -1;
  D3 p, v;
  eval(T, tr, t, p, v, hint);
  double te = t - norm(p - o) / c;
  for (int k = 0; k < 8; ++k) {
    eval(T, tr, te, p, v, hint);
    te = t - norm(p - o) / c;
  }
  if (te < T.t0[tr] || te > T.t1[tr]) return NAN;
  eval(T, tr, te, p, v, hint);
  D3 rv = o - p;
  double vr = dot(v, rv) / norm(rv);
  return (float)(-log1p(-vr / c));
}

// Block-parallel closest approach (scan) + serial golden-section refinement.
__device__ bool block_cpa(const DTracks& T, uint32_t tr, D3 o, double tg,
                          double& t_cpa, double& d, double& speed, float* sd) {
  __shared__ int s_ok;
  __shared__ double s_t, s_d, s_v;
  double a = fmax(tg - 120, T.t0[tr]), b = fmin(tg + 120, T.t1[tr]);
  bool span = b - a >= 2;
  int n = span ? (int)ceil(b - a) + 1 : 0;
  for (int i = threadIdx.x; i < n; i += blockDim.x) {
    int hint = -1;
    D3 p, v;
    eval(T, tr, fmin(b, a + i), p, v, hint);
    D3 r = p - o;
    reinterpret_cast<double*>(sd)[i] = dot(r, r);
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    s_ok = 0;
    if (span) {
      const double* dd = reinterpret_cast<const double*>(sd);
      double best = 1e30;
      int bi = 0;
      for (int i = 0; i < n; ++i)
        if (dd[i] < best) best = dd[i], bi = i;
      if (bi != 0 && bi != n - 1) {
        double bt = fmin(b, a + bi);
        int hint = -1;
        D3 p, v;
        double lo = bt - 1, hi = bt + 1;
        const double g = 0.618033988749895;
        double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
        auto Dd = [&](double t) {
          eval(T, tr, t, p, v, hint);
          D3 r = p - o;
          return dot(r, r);
        };
        double f1 = Dd(x1), f2 = Dd(x2);
        for (int it = 0; it < 30; ++it) {
          if (f1 < f2)
            hi = x2, x2 = x1, f2 = f1, x1 = hi - g * (hi - lo), f1 = Dd(x1);
          else
            lo = x1, x1 = x2, f1 = f2, x2 = lo + g * (hi - lo), f2 = Dd(x2);
        }
        s_t = 0.5 * (lo + hi);
        eval(T, tr, s_t, p, v, hint);
        s_d = norm(p - o);
        s_v = norm(v);
        s_ok = 1;
      }
    }
  }
  __syncthreads();
  t_cpa = s_t, d = s_d, speed = s_v;
  bool ok = s_ok;
  __syncthreads();
  return ok;
}

// ------------------------------------------------- local FP32 trajectories ---
// Per job, the trajectory points around an event are copied to shared memory
// relative to the observer (metres, |x| < ~50 km -> mm precision in FP32) and
// to the event's time reference
// (|t| < ~400 s -> us precision). All inner-loop work is then FP32 on shared
// memory instead of FP64 Hermite evaluation with binary searches in global
// memory (the RTX 3090 runs FP64 at 1/64 rate).
constexpr int kMaxPts = 384;

struct LTraj {  // shared-memory SoA
  float* t;
  float* p;  // xyz interleaved
  float* v;
  int n;
  float t0, t1;  // track coverage relative to the reference time
};

__device__ __forceinline__ void leval(const LTraj& L, float t, float3& p,
                                      float3& v) {
  if (t <= L.t[0]) {
    float dt = t - L.t[0];
    p = make_float3(L.p[0] + L.v[0] * dt, L.p[1] + L.v[1] * dt,
                    L.p[2] + L.v[2] * dt);
    v = make_float3(L.v[0], L.v[1], L.v[2]);
    return;
  }
  int n = L.n;
  if (t >= L.t[n - 1]) {
    const float* q = L.p + 3 * (n - 1);
    const float* w = L.v + 3 * (n - 1);
    float dt = t - L.t[n - 1];
    p = make_float3(q[0] + w[0] * dt, q[1] + w[1] * dt, q[2] + w[2] * dt);
    v = make_float3(w[0], w[1], w[2]);
    return;
  }
  int lo = 0, hi = n;  // upper_bound
  while (lo < hi) {
    int mid = (lo + hi) >> 1;
    if (t < L.t[mid])
      hi = mid;
    else
      lo = mid + 1;
  }
  int i = lo - 1;
  float h = L.t[i + 1] - L.t[i], s = (t - L.t[i]) / h, s2 = s * s, s3 = s2 * s;
  float h00 = 2 * s3 - 3 * s2 + 1, h10 = (s3 - 2 * s2 + s) * h,
        h01 = -2 * s3 + 3 * s2, h11 = (s3 - s2) * h;
  float d00 = (6 * s2 - 6 * s) / h, d10 = 3 * s2 - 4 * s + 1,
        d01 = (-6 * s2 + 6 * s) / h, d11 = 3 * s2 - 2 * s;
  const float *pa = L.p + 3 * i, *va = L.v + 3 * i, *pb = pa + 3, *vb = va + 3;
  p = make_float3(pa[0] * h00 + va[0] * h10 + pb[0] * h01 + vb[0] * h11,
                  pa[1] * h00 + va[1] * h10 + pb[1] * h01 + vb[1] * h11,
                  pa[2] * h00 + va[2] * h10 + pb[2] * h01 + vb[2] * h11);
  v = make_float3(pa[0] * d00 + va[0] * d10 + pb[0] * d01 + vb[0] * d11,
                  pa[1] * d00 + va[1] * d10 + pb[1] * d01 + vb[1] * d11,
                  pa[2] * d00 + va[2] * d10 + pb[2] * d01 + vb[2] * d11);
}
__device__ __forceinline__ float f3dot(float3 a, float3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Loads the window [tref - wpre, tref + wpost] of track `tr` into shared memory
// (relative to o, tref). Returns false if it does not fit (caller falls back to
// the global FP64 path).
__device__ bool load_local(const DTracks& T, uint32_t tr, D3 o, double tref,
                           double wpre, double wpost, LTraj& L, float* st,
                           float* sp, float* sv) {
  __shared__ int s_i0, s_n;
  if (threadIdx.x == 0) {
    uint64_t o0 = T.off[tr];
    int n = (int)(T.off[tr + 1] - o0);
    const double* tt = T.t + o0;
    auto ub = [&](double x) {
      int lo = 0, hi = n;
      while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (x < tt[mid])
          hi = mid;
        else
          lo = mid + 1;
      }
      return lo;
    };
    int i0 = max(0, ub(tref - wpre) - 2), i1 = min(n - 1, ub(tref + wpost) + 1);
    s_i0 = (int)o0 + i0;
    s_n = i1 - i0 + 1;
  }
  __syncthreads();
  int n = s_n, i0 = s_i0;
  if (n > kMaxPts || n < 2) {
    __syncthreads();
    return false;
  }
  for (int i = threadIdx.x; i < n; i += blockDim.x) {
    uint64_t g = (uint64_t)i0 + i;
    st[i] = (float)(T.t[g] - tref);
    sp[3 * i] = (float)(T.p[3 * g] - o.x),
           sp[3 * i + 1] = (float)(T.p[3 * g + 1] - o.y),
           sp[3 * i + 2] = (float)(T.p[3 * g + 2] - o.z);
    sv[3 * i] = (float)T.v[3 * g], sv[3 * i + 1] = (float)T.v[3 * g + 1],
           sv[3 * i + 2] = (float)T.v[3 * g + 2];
  }
  L = {st, sp, sv, n, (float)(T.t0[tr] - tref), (float)(T.t1[tr] - tref)};
  __syncthreads();
  return true;
}

// track_cpa on the local trajectory (observer at origin), same scan +
// golden-section scheme.
__device__ bool local_cpa(const LTraj& L, float tg_lo, float tg_hi,
                          float& t_cpa, float& d, float& speed, float* sd) {
  __shared__ int s_ok;
  __shared__ float s_t, s_d, s_v;
  float a = fmaxf(tg_lo, L.t0), b = fminf(tg_hi, L.t1);
  bool span = b - a >= 2;
  int n = span ? (int)ceilf(b - a) + 1 : 0;
  for (int i = threadIdx.x; i < n; i += blockDim.x) {
    float3 p, v;
    leval(L, fminf(b, a + i), p, v);
    sd[i] = f3dot(p, p);
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    s_ok = 0;
    if (span) {
      float best = 3e38f;
      int bi = 0;
      for (int i = 0; i < n; ++i)
        if (sd[i] < best) best = sd[i], bi = i;
      if (bi != 0 && bi != n - 1) {
        float bt = fminf(b, a + bi), lo = bt - 1, hi = bt + 1;
        const float g = 0.618033988749895f;
        float3 p, v;
        auto D = [&](float t) {
          leval(L, t, p, v);
          return f3dot(p, p);
        };
        float x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo), f1 = D(x1),
              f2 = D(x2);
        for (int it = 0; it < 30; ++it) {
          if (f1 < f2)
            hi = x2, x2 = x1, f2 = f1, x1 = hi - g * (hi - lo), f1 = D(x1);
          else
            lo = x1, x1 = x2, f1 = f2, x2 = lo + g * (hi - lo), f2 = D(x2);
        }
        s_t = 0.5f * (lo + hi);
        leval(L, s_t, p, v);
        s_d = sqrtf(f3dot(p, p));
        s_v = sqrtf(f3dot(v, v));
        s_ok = 1;
      }
    }
  }
  __syncthreads();
  t_cpa = s_t, d = s_d, speed = s_v;
  bool ok = s_ok;
  __syncthreads();
  return ok;
}

// Emission-sampled model (see doppler.cpp model_lng): rows[i] = ln(f/f0)/dln
// for reception times base + i*h (relative to the reference). `tre`/`lge` are
// shared scratch of capacity `cap`.
__device__ void local_model(const LTraj& L, float base, float h, int nfr,
                            float c, double dln, float* rows, float* tre,
                            float* lge, int cap) {
  __shared__ float s_ta;
  __shared__ int s_ne;
  if (threadIdx.x == 0) {
    auto solve = [&](float t) {
      float3 p, v;
      leval(L, t, p, v);
      float te = t - sqrtf(f3dot(p, p)) / c;
      for (int k = 0; k < 8; ++k) {
        leval(L, te, p, v);
        te = t - sqrtf(f3dot(p, p)) / c;
      }
      return te;
    };
    float ta = solve(base) - h, tb = solve(base + (nfr - 1) * h) + h;
    s_ta = ta;
    s_ne = (int)ceilf((tb - ta) / h) + 2;
  }
  __syncthreads();
  const float ta = s_ta;
  const int ne = s_ne;
  if (ne <= cap) {
    for (int j = threadIdx.x; j < ne; j += blockDim.x) {
      float te = ta + j * h;
      float3 p, v;
      leval(L, te, p, v);
      float r = sqrtf(f3dot(p, p));
      tre[j] = te + r / c;
      lge[j] = (te < L.t0 || te > L.t1)
                   ? NAN
                   : -log1pf(f3dot(v, p) / r / c);  // rvec = -p
    }
    __syncthreads();
    for (int i = threadIdx.x; i < nfr; i += blockDim.x) {
      float t = base + i * h;
      int lo = 0, hi = ne - 1;  // last j with tre[j] <= t, clamped to [0, ne-2]
      while (hi - lo > 1) {
        int mid = (lo + hi) >> 1;
        if (tre[mid] <= t)
          lo = mid;
        else
          hi = mid;
      }
      int j = min(lo, ne - 2);
      float w = (t - tre[j]) / (tre[j + 1] - tre[j]);
      float g = (w < 0 || w > 1 || isnan(lge[j]) || isnan(lge[j + 1]))
                    ? NAN
                    : lge[j] + w * (lge[j + 1] - lge[j]);
      rows[i] = (float)((double)g / dln);
    }
  } else {  // window too long for scratch: per-frame fixed-point solve (still
            // FP32 / shared)
    for (int i = threadIdx.x; i < nfr; i += blockDim.x) {
      float t = base + i * h;
      float3 p, v;
      leval(L, t, p, v);
      float te = t - sqrtf(f3dot(p, p)) / c;
      for (int k = 0; k < 8; ++k) {
        leval(L, te, p, v);
        te = t - sqrtf(f3dot(p, p)) / c;
      }
      float g = NAN;
      if (te >= L.t0 && te <= L.t1) {
        leval(L, te, p, v);
        g = -log1pf(f3dot(v, p) / sqrtf(f3dot(p, p)) / c);
      }
      rows[i] = (float)((double)g / dln);
    }
  }
  __syncthreads();
}

__global__ void __launch_bounds__(kBlock, 3)
    k_score(DTracks T, DDem D, const ScoreJob* jobs, int K, int S, double hop,
            double c, double obs, int maxRows, int neCap, ScoreOut* out) {
  // dynamic shared layout (floats): rows[maxRows] | union{ tre,lge[2*neCap] ;
  // combo scores[kCombo] }
  //                                 | local trajectory t[kMaxPts] p[3k] v[3k] |
  //                                 CPA scan 512
  extern __shared__ double sh_d[];
  float* rows = reinterpret_cast<float*>(sh_d);
  float* uni = rows + maxRows;
  float* sc = uni;
  float* tre = uni;
  float* lge = uni + neCap;
  float* st = uni + max(2 * neCap, kCombo);
  float* sp = st + kMaxPts;
  float* sv = sp + 3 * kMaxPts;
  float* sd = sv + 3 * kMaxPts;
  sd += ((sd - reinterpret_cast<float*>(sh_d)) &
         1);  // 8-byte align for the FP64 fallback view
  const ScoreJob J = jobs[blockIdx.x];
  ScoreOut res{};
  D3 o = geo2ecef(J.lat, J.lon, J.elev + obs);
  LTraj L;
  bool local[kMaxEvents];
  PHASE_T0();
  PHASE_CNT(5);
  // pass 1: closest approach for every event (invalid job if any fails)
  for (int k = 0; k < K; ++k) {
    const DEvent& e = c_E[k];
    local[k] = load_local(T, J.track[k], o, J.Tguess[k], e.wpre, e.wpost, L, st,
                          sp, sv);
    if (!local[k]) PHASE_CNT(4);
    PHASE_ADD(0);
    double tc, d, sp_;
    bool ok;
    if (local[k]) {
      float tcf, df, vf;
      ok = local_cpa(L, -120.f, 120.f, tcf, df, vf, sd);
      tc = J.Tguess[k] + tcf, d = df, sp_ = vf;
    } else {
      ok = block_cpa(T, J.track[k], o, J.Tguess[k], tc, d, sp_, sd);
    }
    PHASE_ADD(1);
    if (!ok) {
      if (threadIdx.x == 0) out[blockIdx.x] = res;  // valid = 0
      return;
    }
    res.T[k] = tc + d / c;
    res.d[k] = d;
    res.v[k] = sp_;
  }
  // pass 2: model + alignment. Anchor sweeps shift and f0; others are locked to
  // the anchor's shift.
  float total = 0;
  double dt = 0;
  for (int k = 0; k < K; ++k) {
    const DEvent& e = c_E[k];
    const DEvent& e0 = c_E[0];
    int Sk = k == 0 ? S : 0;
    int nfr = e.fb1 - e.fb0 + 1, n = nfr + 2 * Sk;
    // reception time of row 0, absolute
    double t_first =
        res.T[0] + (e.t0 + (e.fb0 - Sk) * e.hop) - e0.tc - (k ? dt : 0.0);
    if (local[k]) {
      load_local(T, J.track[k], o, J.Tguess[k], e.wpre, e.wpost, L, st, sp, sv);
      PHASE_ADD(0);
      local_model(L, (float)(t_first - J.Tguess[k]), (float)e.hop, n, (float)c,
                  e.dln, rows, tre, lge, neCap);
    } else {
      for (int i = threadIdx.x; i < n; i += blockDim.x) {
        float g = frame_lng(T, J.track[k], o, t_first + i * e.hop, c);
        rows[i] = (float)((double)g / e.dln);
      }
      __syncthreads();
    }
    PHASE_ADD(2);
    int bs, br;
    res.score[k] = block_align(e, rows, nfr, Sk, bs, br, uni, (int)(sd - uni));
    PHASE_ADD(3);
    res.f0[k] = exp(e.lnf0 + br * e.dln);
    if (k == 0) {
      res.shift = bs;
      dt = bs * hop;
    }
    total += res.score[k];
    __syncthreads();
  }
  res.total = total;
  res.valid = 1;
  if (threadIdx.x == 0) out[blockIdx.x] = res;
}

size_t smem_bytes() {
  size_t f = g_max_rows + std::max<size_t>(2 * (size_t)g_ne_cap, kCombo) +
             7 * kMaxPts + 1 + 512;
  return f * sizeof(float);
}

// ========================================================== DSP kernels ===
__global__ void k_frames(const float* x, long n, const float* win, int wl,
                         int nfft, double t0, double hop, double sr, int f0,
                         int nfr, float* buf) {
  long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (idx >= (long)nfr * nfft) return;
  int f = (int)(idx / nfft), i = (int)(idx % nfft);
  float v = 0;
  if (i < wl) {
    long c = llround((t0 + (f0 + f) * hop) * sr);
    long j = c - wl / 2 + i;
    if (j >= 0 && j < n) v = x[j] * win[i];
  }
  buf[idx] = v;
}

__global__ void k_mag(const cufftComplex* X, int nfft, int nfr, int bin0,
                      int nbins, float norm_, float* out) {
  long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
  if (idx >= (long)nfr * nbins) return;
  int f = (int)(idx / nbins), k = (int)(idx % nbins);
  cufftComplex z = X[(size_t)f * (nfft / 2 + 1) + bin0 + k];
  double m = sqrt((double)z.x * z.x + (double)z.y * z.y) * norm_;
  out[idx] = (float)(20.0 * log10(m + 1e-12));
}

// straight_lng in FP32 (the 3090 runs FP64 at 1/64 rate). Only the time
// difference is formed in FP64; for v < c the root expression has no
// cancellation, so FP32 is accurate to ~1e-7 in ln f (grid rows: 1.5e-3).
__device__ __forceinline__ float straight_lng_f(double t_r, float v, float d,
                                                double t0, float c) {
  float T = (float)(t_r - t0);
  float c2 = c * c, a = c2 - v * v;
  float disc = c2 * v * v * T * T + a * d * d;
  float tau = (c2 * T - sqrtf(disc)) / a;
  float x = v * tau, r = sqrtf(x * x + d * d);
  float vr = -v * x / r;
  return -log1pf(-vr / c);
}

__global__ void k_fit(const double* t, const double* lnf, const double* w,
                      int np, const double* vs, int nv, const double* ds,
                      int nd, const double* tcs, int ntc, double c, int side,
                      float* cost, float* lnf0) {
  long idx = blockIdx.x * (long)blockDim.x + threadIdx.x;
  long total = (long)nv * nd * ntc;
  if (idx >= total) return;
  int itc = (int)(idx % ntc), id = (int)((idx / ntc) % nd),
      iv = (int)(idx / ((long)ntc * nd));
  double v = vs[iv], d = ds[id], tc = tcs[itc];
  const double t0 = tc - d / c;
  auto res = [&](int i) {
    return (float)(lnf[i] - (double)straight_lng_f(t[i], (float)v, (float)d, t0,
                                                   (float)c));
  };
  const double cap = 0.02;
  double sw = 0, s = 0;
  for (int i = 0; i < np; ++i) {
    if (side && (t[i] - tc) * side < 0) continue;
    float r = res(i);
    s += w[i] * r;
    sw += w[i];
  }
  if (sw <= 0) {
    cost[idx] = 1e9f;
    lnf0[idx] = 0;
    return;
  }
  double m = s / sw;
  for (int it = 0; it < 2; ++it) {
    double s2 = 0, w2 = 0;
    for (int i = 0; i < np; ++i) {
      if (side && (t[i] - tc) * side < 0) continue;
      float r = res(i);
      if (fabs(r - m) < cap) s2 += w[i] * r, w2 += w[i];
    }
    if (w2 > 0) m = s2 / w2;
  }
  double cs = 0;
  for (int i = 0; i < np; ++i) {
    if (side && (t[i] - tc) * side < 0) continue;
    float r = res(i);
    cs += w[i] * fmin(fabs(r - m), cap);
  }
  cost[idx] = (float)(cs / sw);
  lnf0[idx] = (float)m;
}

template <class T>
struct DBuf {  // grow-only device buffer
  T* p = nullptr;
  size_t n = 0;
  bool ensure(size_t want) {
    if (want <= n) return true;
    if (p) cudaFree(p);
    p = nullptr;
    n = 0;
    if (cudaMalloc(&p, want * sizeof(T)) != cudaSuccess) return false;
    n = want;
    return true;
  }
};

}  // namespace

// ================================================================= API ===
bool compiled() {
  return true;
}
bool enabled() {
  return g_enabled;
}
std::string info() {
  return g_info;
}

bool init(bool want, int dev) {
  g_enabled = false;
  if (!want) {
    g_info = "disabled (--no-gpu)";
    return false;
  }
  int n = 0;
  if (cudaGetDeviceCount(&n) != cudaSuccess || n == 0) {
    g_info = "no CUDA device found";
    return false;
  }
  dev = std::clamp(dev, 0, n - 1);
  cudaDeviceProp pr;
  if (cudaSetDevice(dev) != cudaSuccess ||
      cudaGetDeviceProperties(&pr, dev) != cudaSuccess)
    return false;
  size_t fr = 0, tot = 0;
  cudaMemGetInfo(&fr, &tot);
  char b[256];
  std::snprintf(b, sizeof b, "%s (sm_%d%d, %.1f GB, %.1f GB free)", pr.name,
                pr.major, pr.minor, gb(tot), gb(fr));
  g_info = b;
  g_enabled = true;
  return true;
}

bool stft_db(const float* x, size_t n, double sr, int wl, int nfft, double t0,
             double hop, int nframes, int bin0, int nbins, float* out) {
  const bool fail_value = false;
  if (!g_enabled || nframes <= 0) return false;
  std::vector<float> win(wl);
  double wsum = 0;
  for (int i = 0; i < wl; ++i)
    win[i] = (float)(0.5 - 0.5 * std::cos(2 * kPiD * (i + 0.5) / wl)),
    wsum += win[i];
  float norm_ = (float)(2.0 / wsum);
  float *dx = nullptr, *dwin = nullptr, *dbuf = nullptr, *dout = nullptr;
  cufftComplex* dX = nullptr;
  int batch =
      std::max(1, std::min(nframes, (int)((size_t)256 << 20) / (nfft * 4)));
  CK(cudaMalloc(&dx, n * sizeof(float)));
  CK(cudaMemcpy(dx, x, n * sizeof(float), cudaMemcpyHostToDevice));
  CK(cudaMalloc(&dwin, wl * sizeof(float)));
  CK(cudaMemcpy(dwin, win.data(), wl * sizeof(float), cudaMemcpyHostToDevice));
  CK(cudaMalloc(&dbuf, (size_t)batch * nfft * sizeof(float)));
  CK(cudaMalloc(&dX, (size_t)batch * (nfft / 2 + 1) * sizeof(cufftComplex)));
  CK(cudaMalloc(&dout, (size_t)batch * nbins * sizeof(float)));
  cufftHandle plan;
  int nn[1] = {nfft};
  if (cufftPlanMany(&plan, 1, nn, nullptr, 1, nfft, nullptr, 1, nfft / 2 + 1,
                    CUFFT_R2C, batch) != CUFFT_SUCCESS) {
    cudaFree(dx), cudaFree(dwin), cudaFree(dbuf), cudaFree(dX), cudaFree(dout);
    return false;
  }
  for (int f0 = 0; f0 < nframes; f0 += batch) {
    int nb = std::min(batch, nframes - f0);
    long tot = (long)batch * nfft;
    k_frames<<<(unsigned)((tot + 255) / 256), 256>>>(
        dx, (long)n, dwin, wl, nfft, t0, hop, sr, f0, batch, dbuf);
    cufftExecR2C(plan, dbuf, dX);
    long tm = (long)nb * nbins;
    k_mag<<<(unsigned)((tm + 255) / 256), 256>>>(dX, nfft, nb, bin0, nbins,
                                                 norm_, dout);
    CK(cudaMemcpy(out + (size_t)f0 * nbins, dout, tm * sizeof(float),
                  cudaMemcpyDeviceToHost));
  }
  cufftDestroy(plan);
  cudaFree(dx), cudaFree(dwin), cudaFree(dbuf), cudaFree(dX), cudaFree(dout);
  return true;
}

bool fit_costs(const std::vector<double>& t, const std::vector<double>& lnf,
               const std::vector<double>& w, const std::vector<double>& vs,
               const std::vector<double>& ds, const std::vector<double>& tcs,
               double c, int side, std::vector<float>& cost,
               std::vector<float>& lnf0) {
  const bool fail_value = false;
  if (!g_enabled) return false;
  size_t total = vs.size() * ds.size() * tcs.size();
  if (!total) return false;
  cost.resize(total);
  lnf0.resize(total);
  double *dt, *dl, *dw, *dv, *dd, *dtc;
  float *dc, *df;
  auto up = [&](double*& p, const std::vector<double>& v) {
    return cudaMalloc(&p, v.size() * sizeof(double)) == cudaSuccess &&
           cudaMemcpy(p, v.data(), v.size() * sizeof(double),
                      cudaMemcpyHostToDevice) == cudaSuccess;
  };
  if (!up(dt, t) || !up(dl, lnf) || !up(dw, w) || !up(dv, vs) || !up(dd, ds) ||
      !up(dtc, tcs))
    return false;
  CK(cudaMalloc(&dc, total * sizeof(float)));
  CK(cudaMalloc(&df, total * sizeof(float)));
  k_fit<<<(unsigned)((total + 255) / 256), 256>>>(
      dt, dl, dw, (int)t.size(), dv, (int)vs.size(), dd, (int)ds.size(), dtc,
      (int)tcs.size(), c, side, dc, df);
  CK(cudaGetLastError());
  CK(cudaMemcpy(cost.data(), dc, total * sizeof(float),
                cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(lnf0.data(), df, total * sizeof(float),
                cudaMemcpyDeviceToHost));
  for (void* p : {(void*)dt, (void*)dl, (void*)dw, (void*)dv, (void*)dd,
                  (void*)dtc, (void*)dc, (void*)df})
    cudaFree(p);
  return true;
}

bool upload_tracks(const std::vector<Track>& tracks) {
  const bool fail_value = false;
  if (!g_enabled) return false;
  free_all(g_track_bufs);
  uint64_t N = 0;
  std::vector<uint64_t> off(tracks.size() + 1);
  std::vector<double> t0(tracks.size()), t1(tracks.size());
  for (size_t i = 0; i < tracks.size(); ++i) {
    off[i] = N;
    N += tracks[i].pts.size();
    t0[i] = tracks[i].t0;
    t1[i] = tracks[i].t1;
  }
  off.back() = N;
  double *dt, *dp, *dv, *d0, *d1;
  uint64_t* doff;
  size_t need = N * 7 * sizeof(double) + off.size() * 8 + tracks.size() * 16;
  size_t fr = 0, tot = 0;
  cudaMemGetInfo(&fr, &tot);
  if (need > fr * 0.8) {
    log("gpu: trajectories need %.1f GB but only %.1f GB free on the device - "
        "using CPU",
        gb(need), gb(fr));
    g_enabled = false;
    return false;
  }
  CK(cudaMalloc(&dt, N * sizeof(double)));
  CK(cudaMalloc(&dp, 3 * N * sizeof(double)));
  CK(cudaMalloc(&dv, 3 * N * sizeof(double)));
  CK(cudaMalloc(&doff, off.size() * sizeof(uint64_t)));
  CK(cudaMalloc(&d0, t0.size() * sizeof(double)));
  CK(cudaMalloc(&d1, t1.size() * sizeof(double)));
  g_track_bufs = {dt, dp, dv, doff, d0, d1};
  // stream through a bounded host staging buffer (no full host copy)
  const size_t kStage = 1 << 20;
  std::vector<double> st(kStage), sp(3 * kStage), sv(3 * kStage);
  size_t fill = 0, base = 0;
  auto flush = [&]() -> bool {
    if (!fill) return true;
    bool ok = cudaMemcpy(dt + base, st.data(), fill * 8,
                         cudaMemcpyHostToDevice) == cudaSuccess &&
              cudaMemcpy(dp + 3 * base, sp.data(), fill * 24,
                         cudaMemcpyHostToDevice) == cudaSuccess &&
              cudaMemcpy(dv + 3 * base, sv.data(), fill * 24,
                         cudaMemcpyHostToDevice) == cudaSuccess;
    base += fill;
    fill = 0;
    return ok;
  };
  for (auto& tr : tracks)
    for (auto& q : tr.pts) {
      st[fill] = q.t;
      sp[3 * fill] = q.p.x, sp[3 * fill + 1] = q.p.y, sp[3 * fill + 2] = q.p.z;
      sv[3 * fill] = q.v.x, sv[3 * fill + 1] = q.v.y, sv[3 * fill + 2] = q.v.z;
      if (++fill == kStage && !flush()) return false;
    }
  if (!flush()) return false;
  CK(cudaMemcpy(doff, off.data(), off.size() * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d0, t0.data(), t0.size() * 8, cudaMemcpyHostToDevice));
  CK(cudaMemcpy(d1, t1.data(), t1.size() * 8, cudaMemcpyHostToDevice));
  g_T = {dt, dp, dv, doff, d0, d1, (uint32_t)tracks.size()};
  log("gpu: uploaded %zu trajectories (%llu points, %.2f GB)", tracks.size(),
      (unsigned long long)N, gb(need));
  return true;
}

bool upload_dem(const Dem& dem) {
  const bool fail_value = false;
  if (!g_enabled) return false;
  free_all(g_dem_bufs);
  g_D = DDem{};
  g_D.def = dem.default_elev();
  const auto& tiles = dem.tiles();
  if (tiles.empty()) return true;
  int la0 = 1 << 30, la1 = -(1 << 30), lo0 = 1 << 30, lo1 = -(1 << 30);
  auto dec =
      [](int64_t k, int& la,
         int& lo) {  // inverse of Dem::tkey (la*1000+lo, lo in [-180,180))
        la = (int)std::floor((k + 500.0) / 1000.0);
        lo = (int)(k - (int64_t)la * 1000);
      };
  for (auto& kv : tiles) {
    int la, lo;
    dec(kv.first, la, lo);
    la0 = std::min(la0, la), la1 = std::max(la1, la), lo0 = std::min(lo0, lo),
    lo1 = std::max(lo1, lo);
  }
  int nla = la1 - la0 + 1, nlo = lo1 - lo0 + 1;
  std::vector<int> idx((size_t)nla * nlo, -1), tn;
  std::vector<uint64_t> toff;
  uint64_t tot = 0;
  for (auto& kv : tiles) {
    int la, lo;
    dec(kv.first, la, lo);
    idx[(size_t)(la - la0) * nlo + (lo - lo0)] = (int)tn.size();
    tn.push_back(kv.second.n);
    toff.push_back(tot);
    tot += kv.second.h.size();
  }
  int *didx, *dtn;
  int16_t* ddata;
  uint64_t* dtoff;
  CK(cudaMalloc(&didx, idx.size() * sizeof(int)));
  CK(cudaMalloc(&dtn, tn.size() * sizeof(int)));
  CK(cudaMalloc(&dtoff, toff.size() * 8));
  CK(cudaMalloc(&ddata, tot * sizeof(int16_t)));
  g_dem_bufs = {didx, dtn, dtoff, ddata};
  CK(cudaMemcpy(didx, idx.data(), idx.size() * sizeof(int),
                cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dtn, tn.data(), tn.size() * sizeof(int),
                cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dtoff, toff.data(), toff.size() * 8, cudaMemcpyHostToDevice));
  size_t i = 0;
  for (auto& kv : tiles) {
    CK(cudaMemcpy(ddata + toff[i], kv.second.h.data(), kv.second.h.size() * 2,
                  cudaMemcpyHostToDevice));
    ++i;
  }
  g_D = {didx, ddata, dtoff, dtn, la0, lo0, nla, nlo, dem.default_elev(), true};
  return true;
}

int gen_hits(const HitParams& hp, const std::vector<SampleRange>& ranges_in,
             size_t cap, std::vector<GHit>& out) {
  const int fail_value = -1;
  out.clear();
  if (!g_enabled || !g_T.n) return -1;
  Timer t_all;
  double t_kern = 0;
  // split long ranges (duplicate minima at split points are removed by the
  // dedupe below)
  std::vector<SampleRange> ranges;
  for (auto r : ranges_in)
    for (int32_t k = r.ka; k <= r.kb; k += 65536)
      ranges.push_back({r.track, k, std::min(r.kb, k + 65535), r.a});
  if (ranges.empty()) return 1;
  const double Rh = hp.dmax * 1.03, dlat = hp.step / 111132.954;
  const double mla = Rh / 110000.0 + 0.005;
  const double mlo =
      Rh / (111320.0 * std::cos(std::min(85.0, std::max(std::fabs(hp.la0),
                                                        std::fabs(hp.la1))) *
                                kDegD)) +
      0.005;
  const int maxRows = (int)std::ceil(2 * Rh / hp.step) + 4;
  static DBuf<SampleRange> dR;
  static DBuf<uint64_t> dP;
  static DBuf<Samp> dS;
  static DBuf<Cand> dC;
  static DBuf<GHit> dH;
  static DBuf<unsigned long long> dN;
  const size_t hitCap = std::min<size_t>(cap + 1, (size_t)64 << 20);
  if (!dN.ensure(2) || !dH.ensure(hitCap)) return -1;
  CK(cudaMemset(dN.p + 1, 0, sizeof(unsigned long long)));
  size_t budget = (size_t)1 << 20;  // samples per batch
  size_t candCap = (size_t)40 << 20;
  size_t r0 = 0;
  while (r0 < ranges.size()) {
    // batch of whole ranges up to `budget` samples
    size_t r1 = r0;
    uint64_t N = 0;
    std::vector<uint64_t> prefix;
    while (r1 < ranges.size() &&
           (N == 0 || N + (ranges[r1].kb - ranges[r1].ka + 1) <= budget)) {
      prefix.push_back(N);
      N += ranges[r1].kb - ranges[r1].ka + 1;
      ++r1;
    }
    prefix.push_back(N);
    int nr = (int)(r1 - r0);
    if (!dR.ensure(nr) || !dP.ensure(nr + 1) || !dS.ensure(N) ||
        !dC.ensure(candCap))
      return -1;
    Timer tk;
    CK(cudaMemcpy(dR.p, ranges.data() + r0, nr * sizeof(SampleRange),
                  cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dP.p, prefix.data(), (nr + 1) * 8, cudaMemcpyHostToDevice));
    k_samples<<<(unsigned)((N + 255) / 256), 256>>>(g_T, g_D, dR.p, dP.p, nr, N,
                                                    hp, Rh, mla, mlo, dS.p);
    CK(cudaMemset(dN.p, 0, sizeof(unsigned long long)));
    uint64_t nthr = N * maxRows;
    k_raster<<<(unsigned)((nthr + 255) / 256), 256>>>(
        g_D, dS.p, N, maxRows, hp, Rh, dlat, dC.p, dN.p, candCap);
    unsigned long long nc = 0;
    CK(cudaMemcpy(&nc, dN.p, 8, cudaMemcpyDeviceToHost));
    if (nc >
        candCap) {  // candidate buffer overflow: halve the batch and redo it
      if (budget <= 1024) {
        log("gpu: candidate buffer overflow on a tiny batch");
        return -1;
      }
      budget /= 2;
      continue;
    }
    if (nc) {
      k_refine<<<(unsigned)((nc + 255) / 256), 256>>>(
          g_T, g_D, dS.p, dR.p, dC.p, nc, hp, dlat, dH.p, dN.p + 1, hitCap);
      CK(cudaGetLastError());
    }
    unsigned long long nh = 0;
    CK(cudaMemcpy(&nh, dN.p + 1, 8, cudaMemcpyDeviceToHost));
    t_kern += tk.sec();
    if (nh > cap) return 0;  // over the tile budget -> caller splits the tile
    r0 = r1;
  }
  unsigned long long nh = 0;
  CK(cudaMemcpy(&nh, dN.p + 1, 8, cudaMemcpyDeviceToHost));
  Timer tdl;
  // sort by (node, flight, time) and drop duplicate reports on the device;
  // output is node-sorted, which the time-offset join needs anyway
  size_t w = 0;
  double t_s = 0, t_u = 0, t_r = 0;
  if (nh) {
    thrust::sort(thrust::device, dH.p, dH.p + nh, HitLess());
    cudaDeviceSynchronize();
    t_s = tdl.sec();
    w = thrust::unique(thrust::device, dH.p, dH.p + nh, HitDup()) - dH.p;
    CK(cudaGetLastError());
    t_u = tdl.sec();
  }
  out.resize(w);
  t_r = tdl.sec();
  if (w)
    CK(cudaMemcpy(out.data(), dH.p, w * sizeof(GHit), cudaMemcpyDeviceToHost));
  vlog(
      "gpu gen_hits: %zu ranges, kernels %.3fs, download+dedupe %.3fs, total "
      "%.3fs [sort %.3f unique %.3f resize %.3f copy %.3f, %llu hits]",
      ranges.size(), t_kern, tdl.sec(), t_all.sec(), t_s, t_u - t_s, t_r - t_u,
      tdl.sec() - t_r, nh);
  return 1;
}

bool upload_events(const std::vector<EventImage>& ev, int S, double hop,
                   double c, double obs_h) {
  const bool fail_value = false;
  if (!g_enabled) return false;
  if ((int)ev.size() > kMaxEvents) return false;
  free_all(g_event_bufs);
  g_K = 0;
  const int K_new = (int)ev.size();
  g_S = S;
  g_hop = hop;
  g_c = c;
  g_obs = obs_h;
  g_max_rows = 0;
  for (int k = 0; k < K_new; ++k) {
    const EventImage& e = ev[k];
    if (e.mrow.size() > 8) return false;
    const int nrr = e.rhi >= e.rlo ? (e.rhi - e.rlo) / 3 + 1 : 1,
              Z = ((nrr + 31) / 32) * 32 + 32, L3 = e.nq + 2 * Z;
    const size_t nfr = (size_t)(e.fb1 - e.fb0 + 1);
    size_t n = nfr * e.nb * e.nrH, n3 = nfr * e.nb * 3 * L3;
    float *di, *db;
    CK(cudaMalloc(&di, n * sizeof(float)));
    CK(cudaMalloc(&db, n3 * sizeof(float)));
    g_event_bufs.push_back(di);
    g_event_bufs.push_back(db);
    CK(cudaMemcpy(di, e.himg, n * sizeof(float), cudaMemcpyHostToDevice));
    {
      std::vector<float> pad(n3, 0.f);
      for (size_t fp = 0; fp < nfr * e.nb * 3; ++fp)
        std::copy(e.hblur3 + fp * e.nq, e.hblur3 + (fp + 1) * e.nq,
                  pad.begin() + fp * L3 + Z);
      CK(cudaMemcpy(db, pad.data(), n3 * sizeof(float),
                    cudaMemcpyHostToDevice));
    }
    DEvent d{};
    d.himg = di, d.hblur3 = db, d.nb = e.nb, d.E = e.E, d.nrH = e.nrH,
    d.nq = e.nq, d.L3 = L3, d.Z = Z, d.nf = e.nf, d.nr = e.nr, d.fb0 = e.fb0,
    d.fb1 = e.fb1, d.rlo = e.rlo, d.rhi = e.rhi;
    d.nm = (int)e.mrow.size();
    for (int i = 0; i < 8; ++i)
      d.thr[i] = e.thr[i];
    for (int m = 0; m < d.nm; ++m)
      d.mrow[m] = e.mrow[m];
    d.t0 = e.t0, d.hop = e.hop, d.lnf0 = e.lnf0, d.dln = e.dln, d.tc = e.tc;
    {  // frames relative to this event's CPA; generous margins for node
       // offsets, timing slack, r/c
      double lo = e.t0 + (e.fb0 - (k == 0 ? S : 0)) * e.hop - e.tc,
             hi = e.t0 + (e.fb1 + (k == 0 ? S : 0)) * e.hop - e.tc;
      d.wpre = std::max(130.0, -lo + 100.0) + 6;
      d.wpost = std::max(130.0, hi + 70.0) + 6;
    }
    g_E[k] = d;
    g_max_rows = std::max<size_t>(
        g_max_rows, (size_t)(e.fb1 - e.fb0 + 1) + (k == 0 ? 2 * S : 0));
  }
  CK(cudaMemcpyToSymbol(c_E, g_E, sizeof(DEvent) * kMaxEvents));
  g_K = K_new;
  g_ne_cap = (int)(g_max_rows * 3 / 2 + 64);
  // alignment scratch (offsets + combo scores) reuses the emission + trajectory
  // region; grow if needed
  for (int k = 0; k < K_new; ++k) {
    const int nrow = ev[k].fb1 - ev[k].fb0 + 1 + (k == 0 ? 2 * S : 0);
    const int nrr =
        ev[k].rhi >= ev[k].rlo ? (ev[k].rhi - ev[k].rlo) / 3 + 1 : 1;
    const int need = nrow + std::max(nrr, 128 + nrow);
    const int have = std::max(2 * g_ne_cap, kCombo) + 7 * kMaxPts;
    if (need > have) g_ne_cap += (need - have + 1) / 2;
  }
  size_t smem = smem_bytes();
  if (smem > 48 * 1024) {
    if (cudaFuncSetAttribute(k_score,
                             cudaFuncAttributeMaxDynamicSharedMemorySize,
                             (int)smem) != cudaSuccess) {
      log("gpu: event window too long for shared memory (%zu frames) - scoring "
          "on CPU",
          g_max_rows);
      g_K = 0;
      return false;
    }
  }
  return true;
}

bool score(const std::vector<ScoreJob>& jobs, int K,
           std::vector<ScoreOut>& out) {
  const bool fail_value = false;
  if (!g_enabled || K != g_K || jobs.empty()) return false;
  static DBuf<ScoreJob> dJ;
  static DBuf<ScoreOut> dO;
  out.resize(jobs.size());
  const size_t kChunk = 1 << 16;
  size_t smem = smem_bytes();
  for (size_t b0 = 0; b0 < jobs.size(); b0 += kChunk) {
    size_t n = std::min(kChunk, jobs.size() - b0);
    if (!dJ.ensure(n) || !dO.ensure(n)) return false;
    CK(cudaMemcpy(dJ.p, jobs.data() + b0, n * sizeof(ScoreJob),
                  cudaMemcpyHostToDevice));
    k_score<<<(unsigned)n, kBlock, smem>>>(g_T, g_D, dJ.p, K, g_S, g_hop, g_c,
                                           g_obs, (int)g_max_rows, g_ne_cap,
                                           dO.p);
    CK(cudaGetLastError());
    CK(cudaMemcpy(out.data() + b0, dO.p, n * sizeof(ScoreOut),
                  cudaMemcpyDeviceToHost));
  }
#ifdef PL_PROF_PHASES
  unsigned long long ph[12];
  cudaMemcpyFromSymbol(ph, g_phase, sizeof ph);
  double al = (double)(ph[6] + ph[7] + ph[8] + ph[9]);
  log("  align split: coarse-eval %.1f%%  coarse-serial %.1f%%  refine-eval "
      "%.1f%%  refine-serial %.1f%%",
      100 * ph[6] / al, 100 * ph[7] / al, 100 * ph[8] / al, 100 * ph[9] / al);
  double tot = (double)(ph[0] + ph[1] + ph[2] + ph[3]);
  log("k_score phases (block-cycles): load %.1f%%  cpa %.1f%%  model %.1f%%  "
      "align %.1f%% | jobs %llu, fallback events %llu, "
      "avg %.0f kcycles/job",
      100 * ph[0] / tot, 100 * ph[1] / tot, 100 * ph[2] / tot,
      100 * ph[3] / tot, ph[5], ph[4], tot / std::max(1ULL, ph[5]) / 1e3);
#endif
  return true;
}

}  // namespace pl::gpu
