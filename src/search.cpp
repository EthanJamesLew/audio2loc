// Search pipeline (mirrors the methodology step by step):
//  A. Per event, candidate observer grid nodes: every node whose closest
//  approach (CPA) to a
//     speed-matching, non-jet flight is at the measured slant distance. Each
//     hit stores the *delayed* CPA time (CPA emission time + slant distance /
//     c), using DEM ground elevation.
//  B. Time-offset join: a node must hear event k exactly (tc_k - tc_1) seconds
//  after event 1. C. Silence filter: the node must NOT have heard any other
//  aircraft in the gaps. D. Full-geometry Doppler modelling on a finer grid:
//  for the actual ADS-B trajectories, compute
//     the curve each node would have heard (retarded-time radial velocity), and
//     score it by summing spectrogram brightness under it. Event 1 sweeps
//     time/f0 to absorb measurement error; the other events are LOCKED to event
//     1 by flight geometry and only sweep f0.
//  E. Group by flight combination, keep the best node per group, rank; dense
//  heatmap for the winner.
#include "search.h"

#include <array>
#include <future>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

#include "export.h"
#include "gpu.h"
#include "util.h"

namespace fs = std::filesystem;

namespace pl {

// ------------------------------------------------------- event prep ---
namespace {
void build_harmonic_tables(Event& e, const std::vector<float>& mrow);
}

std::vector<Event> prepare_events(const Audio& a,
                                  const std::vector<EventSpec>& specs,
                                  const SearchParams& P) {
  std::vector<Event> out;
  int S = (int)std::ceil(P.dt_sweep / P.hop);
  for (size_t i = 0; i < specs.size(); ++i) {
    const EventSpec& s = specs[i];
    Event e;
    e.spec = s;
    double mmin = *std::min_element(s.mult.begin(), s.mult.end());
    double mmax = *std::max_element(s.mult.begin(), s.mult.end());
    double pad = S * P.hop + 20;
    double flo = std::min(s.fmin, s.fmin * mmin) * 0.8,
           fhi = std::max(s.fmax, s.fmax * mmax) * 1.2;
    Spectrogram sp =
        compute_spectrogram(a, P.win, P.hop, s.t0 - pad, s.t1 + pad, flo, fhi);
    if (sp.nframes < 10)
      die("event %zu: box [%g,%g] s is outside the audio (%.1f s)", i + 1, s.t0,
          s.t1, a.duration());
    whiten(sp);
    Ridge full;  // continuity-constrained path with no brightness threshold
                 // (used for ridge scoring)
    if (!s.trace.empty()) {
      std::ifstream tf(s.trace);
      if (!tf) die("cannot read trace %s", s.trace.c_str());
      std::string line;
      while (std::getline(tf, line)) {
        auto c = split(line, line.find(',') != std::string::npos ? ',' : ' ');
        double t, f;
        if (c.size() >= 2 && parse_double(c[0], t) && parse_double(c[1], f) &&
            t >= s.t0 && t <= s.t1)
          e.ridge.t.push_back(t), e.ridge.f.push_back(f),
              e.ridge.w.push_back(5.0);
      }
      log("event %zu: %zu traced points from %s", i + 1, e.ridge.t.size(),
          s.trace.c_str());
      full = e.ridge;
    } else {
      e.ridge = extract_ridge(sp, s.t0, s.t1, s.fmin, s.fmax);
      full = extract_ridge(sp, s.t0, s.t1, s.fmin, s.fmax, 25.0, 0.0f);
    }
    e.fit = fit_straight(e.ridge, P.c, s.side);
    auto pick = [](double over, double meas) {
      return std::isnan(over) ? meas : over;
    };
    e.vmin = pick(s.vmin, e.fit.vmin * (1 - P.v_margin));
    e.vmax = pick(s.vmax, e.fit.vmax * (1 + P.v_margin));
    e.dmin = pick(s.dmin, e.fit.dmin * (1 - P.d_margin));
    e.dmax = pick(s.dmax, e.fit.dmax * (1 + P.d_margin));
    e.tc = pick(s.tc, e.fit.tc);
    e.tc_tol = pick(
        s.tc_tol,
        std::max(e.fit.tcmax - e.fit.tc, e.fit.tc - e.fit.tcmin) + P.tc_margin);
    e.f0min = pick(s.f0min, e.fit.f0min * (1 - P.f0_margin));
    e.f0max = pick(s.f0max, e.fit.f0max * (1 + P.f0_margin));
    e.white = sp;
    e.img = to_log_image(sp, flo, fhi, P.dln);
    if (s.ridge_score) {
      // Faint-curve mode: score by closeness to the traced ridge instead of raw
      // brightness.
      std::fill(e.img.v.begin(), e.img.v.end(), 0.f);
      double sig = 0.003 / P.dln;
      for (size_t k = 0; k < full.t.size(); ++k) {
        int f = e.img.frame_of(full.t[k]);
        double r0 = e.img.row_of(full.f[k]);
        for (int r = (int)(r0 - 4 * sig); r <= (int)(r0 + 4 * sig) + 1; ++r)
          if (f >= 0 && f < e.img.nf && r >= 0 && r < e.img.nr)
            e.img.v[(size_t)f * e.img.nr + r] =
                (float)(8.0 * std::exp(-(r - r0) * (r - r0) / (2 * sig * sig)));
      }
      e.spec.mult = {1.0};
    }
    e.blur = e.img.blurred(1, 2);
    e.fb0 = std::clamp(e.img.frame_of(s.t0), 0, e.img.nf - 1);
    e.fb1 = std::clamp(e.img.frame_of(s.t1), 0, e.img.nf - 1);
    {
      std::vector<float> mrow;
      for (double k : e.spec.mult)
        mrow.push_back((float)(std::log(k) / e.img.dln));
      build_harmonic_tables(e, mrow);
    }
    out.push_back(std::move(e));
  }
  return out;
}

void print_events(const std::vector<Event>& ev) {
  for (size_t i = 0; i < ev.size(); ++i) {
    const Event& e = ev[i];
    const DopplerFit& F = e.fit;
    std::printf(
        "event %zu  box %.1f-%.1f s, %.0f-%.0f Hz, %d ridge pts\n"
        "  straight-line fit : v=%.1f m/s (%.0f kt)  d_cpa=%.0f m  tc=%.2f s  "
        "f0=%.2f Hz  cost=%.4f\n"
        "  fit uncertainty   : v %.1f-%.1f  d %.0f-%.0f  tc %.2f-%.2f  f0 "
        "%.2f-%.2f\n"
        "  search ranges     : v %.1f-%.1f m/s  d %.0f-%.0f m  tc %.2f+-%.2f s "
        " f0 %.2f-%.2f Hz\n",
        i + 1, e.spec.t0, e.spec.t1, e.spec.fmin, e.spec.fmax, F.npts, F.v,
        F.v / kKt, F.d, F.tc, F.f0, F.cost, F.vmin, F.vmax, F.dmin, F.dmax,
        F.tcmin, F.tcmax, F.f0min, F.f0max, e.vmin, e.vmax, e.dmin, e.dmax,
        e.tc, e.tc_tol, e.f0min, e.f0max);
  }
  if (ev.size() > 1)
    for (size_t i = 1; i < ev.size(); ++i)
      std::printf("offset event %zu - event 1: %.2f s\n", i + 1,
                  ev[i].tc - ev[0].tc);
}

namespace {

struct Hit {
  uint64_t key;
  uint32_t track;
  float d, v;
  double T;  // delayed CPA (reception) absolute time
};

// ---------------------------------------------------- stage A: hits ---
struct Tile {
  double la0, la1, lo0, lo1;
  bool owns(double la, double lo) const {
    return la >= la0 && la < la1 && lo >= lo0 && lo < lo1;
  }
};

// Time intervals during which a track can be within `mla`/`mlo` degrees of the
// tile.
static void tile_intervals(const Track& tr, const Tile& T, double mla,
                           double mlo,
                           std::vector<std::pair<double, double>>& iv) {
  iv.clear();
  const auto& p = tr.pts;
  const double pad = 0.02;  // Hermite bulge beyond the chord
  for (size_t k = 0; k + 1 < p.size(); ++k) {
    double a0 = std::min(p[k].lat, p[k + 1].lat) - pad,
           a1 = std::max(p[k].lat, p[k + 1].lat) + pad;
    double o0 = std::min(p[k].lon, p[k + 1].lon) - pad,
           o1 = std::max(p[k].lon, p[k + 1].lon) + pad;
    if (a1 < T.la0 - mla || a0 > T.la1 + mla || o1 < T.lo0 - mlo ||
        o0 > T.lo1 + mlo)
      continue;
    if (!iv.empty() && p[k].t <= iv.back().second)
      iv.back().second = p[k + 1].t;
    else
      iv.push_back({p[k].t, p[k + 1].t});
  }
}

// CPA hits for grid nodes OWNED by `tile`. Per-thread scratch is bounded by
// processing each track in chunks of `chunk` samples. Returns false (and stops
// early) once more than `cap` hits are produced: the caller then splits the
// tile, which is what keeps memory under the budget.
bool gen_hits(const Event& e, const std::vector<Track>& tracks,
              const std::vector<uint32_t>& subset, const GeoGrid& G,
              const Dem& dem, const SearchParams& P, double Tlo, double Thi,
              const Tile& tile, size_t cap, int chunk, std::vector<Hit>& out) {
  static_assert(sizeof(Hit) == sizeof(gpu::GHit),
                "Hit / GHit layout must match");
  if (gpu::enabled()) {
    // Same sample set as the CPU loop below, described as (track, k-range) for
    // the kernels.
    const double c = P.c, dt = P.sample_dt, Rh = e.dmax * 1.03;
    const double mla = Rh / 110000.0 + 0.005;
    const double mlo =
        Rh /
            (111320.0 * std::cos(std::min(85.0, std::max(std::fabs(tile.la0),
                                                         std::fabs(tile.la1))) *
                                 kDeg)) +
        0.005;
    std::vector<gpu::SampleRange> ranges;
    std::vector<std::pair<double, double>> iv;
    for (uint32_t ti : subset) {
      const Track& tr = tracks[ti];
      double a = std::max(tr.t0, Tlo - e.dmax / c - 5),
             b = std::min(tr.t1, Thi + 5);
      if (b - a < 2 * dt) continue;
      tile_intervals(tr, tile, mla, mlo, iv);
      int kmax = (int)std::floor((b - a) / dt + 1e-9);
      for (auto [ia, ib] : iv) {
        ia = std::max(ia, a);
        ib = std::min(ib, b);
        if (ib <= ia) continue;
        int ka = std::max(0, (int)std::floor((ia - a) / dt)),
            kb = std::min(kmax, (int)std::ceil((ib - a) / dt));
        if (kb >= ka) ranges.push_back({ti, ka, kb, a});
      }
    }
    vlog("gen_hits host prep: %zu ranges", ranges.size());
    gpu::HitParams hp{e.vmin,   e.vmax,   e.dmin,   e.dmax,       Tlo,
                      Thi,      c,        dt,       P.obs_height, G.step(),
                      tile.la0, tile.la1, tile.lo0, tile.lo1};
    std::vector<gpu::GHit> gh;
    int rc = gpu::gen_hits(hp, ranges, cap, gh);
    if (rc >= 0) {
      out.resize(gh.size());
      if (!gh.empty())
        std::memcpy(out.data(), gh.data(), gh.size() * sizeof(Hit));
      return rc == 1;
    }
    // rc < 0: GPU failed -> CPU path below
  }
  int nt = num_threads();
  std::vector<std::vector<Hit>> part(nt);
  struct TH {
    uint64_t key;
    int k;
    float d;
  };
  std::vector<std::vector<TH>> scratch(nt);
  std::vector<std::vector<std::pair<double, double>>> ivs(nt);
  std::atomic<size_t> total{0};
  std::atomic<bool> over{false};
  const double c = P.c, dt = P.sample_dt, Rh = e.dmax * 1.03;
  const double mla = Rh / 110000.0 + 0.005;
  const double mlo =
      Rh / (111320.0 * std::cos(std::min(85.0, std::max(std::fabs(tile.la0),
                                                        std::fabs(tile.la1))) *
                                kDeg)) +
      0.005;
  parallel_for(subset.size(), [&](size_t ii, int tid) {
    if (over.load(std::memory_order_relaxed)) return;
    const uint32_t ti = subset[ii];
    const Track& tr = tracks[ti];
    double a = std::max(tr.t0, Tlo - e.dmax / c - 5),
           b = std::min(tr.t1, Thi + 5);
    if (b - a < 2 * dt) return;
    auto& iv = ivs[tid];
    tile_intervals(tr, tile, mla, mlo, iv);
    auto& th = scratch[tid];
    int hint = -1;
    Vec3 p, v;
    for (auto [ia, ib] : iv) {
      ia = std::max(ia, a);
      ib = std::min(ib, b);
      if (ib <= ia) continue;
      // sample indices are global (relative to `a`) so runs/minima are
      // consistent across intervals
      int ka = std::max(0, (int)std::floor((ia - a) / dt)),
          kb = (int)std::ceil((ib - a) / dt);
      for (int k0 = ka; k0 <= kb; k0 += chunk) {
        int k1 = std::min(
            kb + 1,
            k0 +
                chunk);  // owned [k0, k1); rasterize one extra sample each side
        th.clear();
        for (int k = std::max(ka, k0 - 1); k <= std::min(kb, k1); ++k) {
          double t = a + k * dt;
          if (t > b) break;
          tr.eval(t, p, v, &hint);
          double sp = v.norm();
          if (sp < e.vmin * 0.97 || sp > e.vmax * 1.03) continue;
          double lat, lon, h;
          ecef_to_geodetic(p, lat, lon, h);
          if (lat < tile.la0 - mla || lat > tile.la1 + mla ||
              lon < tile.lo0 - mlo || lon > tile.lo1 + mlo)
            continue;
          if (h - dem.elev(lat, lon) - 600 > Rh) continue;
          double mlat = m_per_deg_lat(lat), mlon = m_per_deg_lon(lat);
          G.for_each_within(
              lat, lon, Rh, [&](uint64_t key, double la, double lo) {
                if (!tile.owns(la, lo)) return;
                double dz = h - (dem.elev(la, lo) + P.obs_height);
                double dx = (lo - lon) * mlon, dy = (la - lat) * mlat;
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 <= Rh * Rh) th.push_back({key, k, (float)std::sqrt(d2)});
              });
        }
        if (th.empty()) continue;
        std::sort(th.begin(), th.end(), [](const TH& x, const TH& y) {
          return x.key != y.key ? x.key < y.key : x.k < y.k;
        });
        size_t i = 0;
        while (i < th.size()) {
          size_t j = i;
          while (j < th.size() && th[j].key == th[i].key)
            ++j;
          double lastT = -1e300;
          for (size_t q = i; q < j; ++q) {
            if (th[q].k < k0 || th[q].k >= k1)
              continue;  // neighbour samples only inform the test
            bool lmin = (q == i || th[q - 1].k != th[q].k - 1 ||
                         th[q - 1].d >= th[q].d) &&
                        (q + 1 == j || th[q + 1].k != th[q].k + 1 ||
                         th[q + 1].d > th[q].d);
            if (!lmin) continue;
            double la, lo;
            G.center(th[q].key, la, lo);
            Vec3 o = geodetic_to_ecef(la, lo, dem.elev(la, lo) + P.obs_height);
            double tk = a + th[q].k * dt, tcpa, d, sp;
            if (!track_cpa(tr, o, tk - 2 * dt - 1, tk + 2 * dt + 1, tcpa, d,
                           sp))
              continue;
            if (d < e.dmin || d > e.dmax || sp < e.vmin || sp > e.vmax)
              continue;
            double T = tcpa + d / c;
            if (T < Tlo || T > Thi || std::fabs(T - lastT) < 2) continue;
            lastT = T;
            part[tid].push_back({th[q].key, ti, (float)d, (float)sp, T});
            if (total.fetch_add(1, std::memory_order_relaxed) + 1 > cap) {
              over = true;
              return;
            }
          }
          i = j;
        }
      }
    }
    if (th.capacity() > (size_t)chunk * 4096)
      std::vector<TH>().swap(th);  // never hoard a huge buffer
  });
  out.clear();
  if (over) return false;
  out.reserve(total.load());
  for (auto& v : part) {
    out.insert(out.end(), v.begin(), v.end());
    std::vector<Hit>().swap(v);
  }
  return true;
}

struct Candidate {
  uint64_t key;
  std::vector<uint32_t> track;
  std::vector<double> T;
};

// ----------------------------------------------- stage C: silence ---
// Set of whole minutes (epoch / 60) as a flat bitmap over [lo, lo + n).
struct MinuteSet {
  int64_t lo = 0;
  std::vector<char> b;
  void build(int64_t mlo, int64_t mhi) {
    lo = mlo, b.assign((size_t)std::max<int64_t>(0, mhi - mlo + 1), 0);
  }
  void insert(int64_t m) { b[(size_t)(m - lo)] = 1; }
  size_t count(int64_t m) const {
    return m >= lo && m < lo + (int64_t)b.size() && b[(size_t)(m - lo)];
  }
};

class SilenceIndex {
 public:
  struct Seg {
    double ta, tb;
    Vec3 A, B;
    uint32_t track;
    bool jet;
  };
  SilenceIndex(const std::vector<Track>& tracks,
               const std::vector<uint32_t>& subset, const MinuteSet& minutes,
               double step, double la0 = -90, double la1 = 90,
               double lo0 = -360, double lo1 = 360) {
    // only cells a query from inside [la0,la1]x[lo0,lo1] can visit (dy +-1, dx
    // +-2)
    const int cy0 = cell(la0) - 1, cy1 = cell(la1) + 1, cx0 = cell(lo0) - 2,
              cx1 = cell(lo1) + 2;
    // Segments are built and bucket-keyed in parallel, then sorted by bucket so
    // every bucket is one contiguous range of segs_ (CSR); the hash map only
    // stores [begin, end).
    std::vector<std::vector<std::pair<uint64_t, Seg>>> part(subset.size());
    parallel_for(
        subset.size(),
        [&](size_t ii, int) {
          const uint32_t ti = subset[ii];
          const Track& tr = tracks[ti];
          int hint = -1;
          Vec3 pa, va, pb, vb;
          for (double t = tr.t0; t < tr.t1; t += step) {
            double tb = std::min(tr.t1, t + step);
            int64_t m = (int64_t)std::floor(t / 60);
            if (!minutes.count(m) && !minutes.count(m + 1) &&
                !minutes.count(m - 1))
              continue;
            tr.eval(t, pa, va, &hint);
            tr.eval(tb, pb, vb, &hint);
            double lat, lon, h;
            ecef_to_geodetic((pa + pb) * 0.5, lat, lon, h);
            const int cy = cell(lat), cx = cell(lon);
            if (cy < cy0 || cy > cy1 || cx < cx0 || cx > cx1) continue;
            part[ii].push_back(
                {bkey(m, cell(lat), cell(lon)),
                 {t, tb, pa, pb, (uint32_t)ti, tr.cls == ACClass::Jet}});
          }
        },
        8);
    std::vector<std::pair<uint64_t, Seg>> all;
    size_t tot = 0;
    for (auto& v : part)
      tot += v.size();
    all.reserve(tot);
    for (auto& v : part)
      all.insert(all.end(), v.begin(), v.end()),
          std::vector<std::pair<uint64_t, Seg>>().swap(v);
    std::stable_sort(all.begin(), all.end(), [](const auto& x, const auto& y) {
      return x.first < y.first;
    });
    segs_.reserve(all.size());
    for (size_t i = 0; i < all.size(); ++i) {
      if (!i || all[i].first != all[i - 1].first)
        buckets_[all[i].first] = {(uint32_t)i, (uint32_t)i};
      buckets_[all[i].first].second = (uint32_t)i + 1;
      segs_.push_back(all[i].second);
    }
    vlog("silence index: %zu segments in %zu buckets", segs_.size(),
         buckets_.size());
  }
  // True if any aircraft (other than excluded event passes) is audible within
  // (Ta, Tb).
  template <class Excl>
  bool heard(const Vec3& o, double lat, double lon, double Ta, double Tb,
             double Rp, double Rj, double c, Excl&& excluded) const {
    double R = std::max(Rp, Rj);
    int64_t m0 = (int64_t)std::floor((Ta - R / c - 10) / 60),
            m1 = (int64_t)std::floor(Tb / 60);
    int cy = cell(lat), cx = cell(lon);
    for (int64_t m = m0; m <= m1; ++m)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -2; dx <= 2; ++dx) {
          auto it = buckets_.find(bkey(m, cy + dy, cx + dx));
          if (it == buckets_.end()) continue;
          for (uint32_t si = it->second.first; si < it->second.second; ++si) {
            const Seg& s = segs_[si];
            Vec3 ab = s.B - s.A;
            double L2 = ab.norm2();
            double u =
                L2 > 0 ? std::clamp((o - s.A).dot(ab) / L2, 0.0, 1.0) : 0;
            double d = (o - (s.A + ab * u)).norm();
            if (d > (s.jet ? Rj : Rp)) continue;
            double trx = s.ta + u * (s.tb - s.ta) + d / c;
            if (trx <= Ta || trx >= Tb) continue;
            if (excluded(s.track, trx)) continue;
            return true;
          }
        }
    return false;
  }

 private:
  static int cell(double deg) { return (int)std::floor(deg / 0.05); }
  static uint64_t bkey(int64_t m, int cy, int cx) {
    return ((uint64_t)m << 25) ^ ((uint64_t)(cy + 2048) << 13) ^
           (uint64_t)(cx + 4096);
  }
  std::vector<Seg> segs_;
  std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> buckets_;
};

// --------------------------------------------- stage D: scoring ---
// rows[i] (i in [0, nfr+2S)) is the model log-frequency row offset ln(g)/dln
// for frame fb0-S+i. Shift s in [-S,S] moves the curve by s frames; r is the f0
// row. Returns mean brightness.
void build_harmonic_tables(Event& e, const std::vector<float>& mrow) {
  auto& H = e.H;
  const int nm = (int)mrow.size(), nr = e.img.nr;
  std::vector<float> th;
  for (float m : mrow) {
    float t = 0.5f - (m - std::floor(m));
    if (t < 0) t += 1.f;
    if (t > 0.f && t < 1.f) th.push_back(t);
  }
  std::sort(th.begin(), th.end());
  th.erase(std::unique(th.begin(), th.end()), th.end());
  H.nb = (int)th.size() + 1;
  for (size_t i = 0; i < th.size(); ++i)
    H.thr[i] = th[i];
  std::vector<std::vector<int>> off(H.nb, std::vector<int>(nm));
  int omax = 0, omin = 0;
  for (int b = 0; b < H.nb; ++b) {
    float lo = b ? th[b - 1] : 0.f, hi = b < H.nb - 1 ? th[b] : 1.f,
          phi = 0.5f * (lo + hi);
    for (int m = 0; m < nm; ++m) {
      off[b][m] = (int)std::lrint(phi + mrow[m]);
      omax = std::max(omax, off[b][m]), omin = std::min(omin, off[b][m]);
    }
  }
  H.E = omax;  // rows below 0 whose upper harmonics still land in the image
  H.nrH = nr + omax - omin;
  H.nq = (H.nrH + 2) / 3;
  H.nfr = e.fb1 - e.fb0 + 1;
  H.himg.assign((size_t)H.nfr * H.nb * H.nrH, 0.f);
  H.hblur3.assign((size_t)H.nfr * H.nb * 3 * H.nq, 0.f);
  parallel_for((size_t)H.nfr, [&](size_t j, int) {
    const int f = e.fb0 + (int)j;
    for (int b = 0; b < H.nb; ++b) {
      float* Hi = &H.himg[(j * H.nb + b) * H.nrH];
      float* Hb = &H.hblur3[(j * H.nb + b) * 3 * H.nq];
      for (int rr = 0; rr < H.nrH; ++rr) {
        float si = 0, sb = 0;
        for (int m = 0; m < nm; ++m) {
          int r = rr - H.E + off[b][m];
          si += e.img.at(f, r), sb += e.blur.at(f, r);
        }
        Hi[rr] = si;
        Hb[(rr % 3) * H.nq + rr / 3] = sb;
      }
    }
  });
}

// Alignment sweep on the harmonic tables (see Event::Harm). For a fixed shift s
// every f0 offset r shares the per-frame (variant, floor row), so the sweep
// over r is a contiguous vector add (row-phase-split table for the step-3
// coarse grid, 7-wide windows for the refine). Same search order and tie rules
// as the original per-combo evaluation.
float best_alignment(const Event& e, const float* rows, int nfr, int S, int fb0,
                     int r_lo, int r_hi, const std::vector<float>& mrow,
                     int& bs, int& br) {
  (void)fb0;
  struct C {
    float sc;
    int s, r;
  };
  const auto& H = e.H;
  C top[3] = {{-1, 0, r_lo}, {-1, 0, r_lo}, {-1, 0, r_lo}};
  const int ss = S >= 2 ? 2 : 1, nm = (int)mrow.size(), nq = H.nq, nb = H.nb,
            nrH = H.nrH;
  const int nrr = r_hi >= r_lo ? (r_hi - r_lo) / 3 + 1 : 0;
  // per model row: table variant and floor row (shared by every shift / r)
  thread_local std::vector<int> vb, vf;
  const int nrow = nfr + 2 * S;
  vb.resize(nrow), vf.resize(nrow);
  for (int i = 0; i < nrow; ++i)
    vb[i] = std::isnan(rows[i]) ? -1 : H.variant(rows[i], vf[i]);
  thread_local std::vector<float> acc;
  acc.resize(nrr);
  const float norm_c = (float)(((nfr + 2) / 3) * nm);
  for (int s = -S; s <= S; s += ss) {
    std::fill(acc.begin(), acc.end(), 0.f);
    for (int j = 0; j < nfr; j += 3) {
      const int i = j + S - s, b = vb[i];
      if (b < 0) continue;
      const float* P = H.hblur3.data() + ((size_t)j * nb + b) * 3 * nq;
      const int row0 = vf[i] + r_lo + H.E;
      const int p = ((row0 % 3) + 3) % 3, q0 = (row0 - p) / 3;
      const float* Q = P + (size_t)p * nq + q0;
      const int ka = std::max(0, -q0), kb = std::min(nrr, nq - q0);
      float* __restrict A = acc.data();
      for (int k = ka; k < kb; ++k)
        A[k] += Q[k];
    }
    for (int k = 0; k < nrr; ++k) {
      float sc = acc[k] / norm_c;
      if (sc > top[2].sc) {
        top[2] = {sc, s, r_lo + 3 * k};
        std::sort(top, top + 3,
                  [](const C& a, const C& b) { return a.sc > b.sc; });
      }
    }
  }
  float best = -1;
  bs = 0;
  br = r_lo;
  const float norm_f = (float)(nfr * nm);
  for (auto& t : top) {
    if (t.sc < 0) continue;
    for (int s = std::max(-S, t.s - ss); s <= std::min(S, t.s + ss); ++s) {
      int ra = std::max(r_lo, t.r - 3), rb = std::min(r_hi, t.r + 3),
          w = rb - ra + 1;
      float a7[7] = {0, 0, 0, 0, 0, 0, 0};
      for (int j = 0; j < nfr; ++j) {
        const int i = j + S - s, b = vb[i];
        if (b < 0) continue;
        const float* P = H.himg.data() + ((size_t)j * nb + b) * nrH;
        const int row0 = vf[i] + ra + H.E;
        for (int q = 0; q < w; ++q)
          if ((unsigned)(row0 + q) < (unsigned)nrH) a7[q] += P[row0 + q];
      }
      for (int q = 0; q < w; ++q) {
        float sc = a7[q] / norm_f;
        if (sc > best) best = sc, bs = s, br = ra + q;
      }
    }
  }
  return best;
}

struct Scorer {
  const std::vector<Track>& tracks;
  const std::vector<Event>& ev;
  const Dem& dem;
  const SearchParams& P;
  int S;
  std::vector<std::vector<float>> mrow;
  std::vector<int> rlo, rhi;
  Scorer(const std::vector<Track>& t, const std::vector<Event>& e, const Dem& d,
         const SearchParams& p)
      : tracks(t), ev(e), dem(d), P(p) {
    S = (int)std::ceil(P.dt_sweep / P.hop);
    for (auto& x : ev) {
      std::vector<float> m;
      for (double k : x.spec.mult)
        m.push_back((float)(std::log(k) / x.img.dln));
      mrow.push_back(m);
      rlo.push_back((int)std::floor(x.img.row_of(x.f0min)));
      rhi.push_back((int)std::ceil(x.img.row_of(x.f0max)));
    }
  }
  // Scores one observer location against one flight combination. Tguess:
  // approximate delayed CPA times.
  bool score(double lat, double lon, const std::vector<uint32_t>& trk,
             const std::vector<double>& Tguess, Scored& out,
             std::vector<double>& tbuf, std::vector<float>& gbuf) const {
    size_t K = ev.size();
    out.lat = lat;
    out.lon = lon;
    out.elev = dem.elev(lat, lon);
    out.track = trk;
    out.T.assign(K, 0);
    out.d.assign(K, 0);
    out.v.assign(K, 0);
    out.f0.assign(K, 0);
    out.score.assign(K, 0);
    Vec3 o = geodetic_to_ecef(lat, lon, out.elev + P.obs_height);
    for (size_t k = 0; k < K; ++k) {
      double tcpa, d, sp;
      if (!track_cpa(tracks[trk[k]], o, Tguess[k] - 120, Tguess[k] + 120, tcpa,
                     d, sp))
        return false;
      out.T[k] = tcpa + d / P.c;
      out.d[k] = d;
      out.v[k] = sp;
    }
    // anchor event: sweep shift and f0
    const Event& e0 = ev[0];
    int nfr0 = e0.fb1 - e0.fb0 + 1, n0 = nfr0 + 2 * S;
    tbuf.resize(n0);
    gbuf.resize(n0);
    for (int i = 0; i < n0; ++i)
      tbuf[i] = out.T[0] + e0.img.time(e0.fb0 - S + i) - e0.tc;
    model_lng(tracks[trk[0]], o, tbuf.data(), n0, P.c, gbuf.data());
    for (auto& g : gbuf)
      g = (float)(g / e0.img.dln);
    int bs, br;
    out.score[0] = best_alignment(e0, gbuf.data(), nfr0, S, e0.fb0, rlo[0],
                                  rhi[0], mrow[0], bs, br);
    out.shift = bs;
    out.f0[0] = e0.img.freq(br);
    double dt = bs * P.hop;
    out.total = out.score[0];
    // locked events: placement fixed by geometry relative to event 1; sweep f0
    // only
    for (size_t k = 1; k < K; ++k) {
      const Event& e = ev[k];
      int nfr = e.fb1 - e.fb0 + 1;
      tbuf.resize(nfr);
      gbuf.resize(nfr);
      for (int i = 0; i < nfr; ++i)
        tbuf[i] = out.T[0] + e.img.time(e.fb0 + i) - e0.tc - dt;
      model_lng(tracks[trk[k]], o, tbuf.data(), nfr, P.c, gbuf.data());
      for (auto& g : gbuf)
        g = (float)(g / e.img.dln);
      int s2, r2;
      out.score[k] = best_alignment(e, gbuf.data(), nfr, 0, e.fb0, rlo[k],
                                    rhi[k], mrow[k], s2, r2);
      out.f0[k] = e.img.freq(r2);
      out.total += out.score[k];
    }
    return true;
  }
  // Batch scoring: GPU when available, else CPU threads. `lat/lon/trk/Tg` per
  // job.
  void score_many(const std::vector<std::pair<double, double>>& ll,
                  const std::vector<const std::vector<uint32_t>*>& trk,
                  const std::vector<const std::vector<double>*>& Tg,
                  std::vector<Scored>& out, std::vector<char>& valid,
                  std::vector<std::vector<double>>& tb,
                  std::vector<std::vector<float>>& gb) const {
    size_t n = ll.size(), K = ev.size();
    out.assign(n, Scored{});
    valid.assign(n, 0);
    if (gpu::enabled() && K <= (size_t)gpu::kMaxEvents) {
      std::vector<gpu::ScoreJob> jobs(n);
      parallel_for(
          n,
          [&](size_t i, int) {
            auto& j = jobs[i];
            j.lat = ll[i].first, j.lon = ll[i].second,
            j.elev = dem.elev(j.lat, j.lon);
            for (size_t k = 0; k < K; ++k)
              j.track[k] = (*trk[i])[k], j.Tguess[k] = (*Tg[i])[k];
          },
          4096);
      std::vector<gpu::ScoreOut> go;
      if (gpu::score(jobs, (int)K, go)) {
        parallel_for(
            n,
            [&](size_t i, int) {
              const auto& g = go[i];
              if (!g.valid) return;
              Scored& s = out[i];
              s.lat = jobs[i].lat, s.lon = jobs[i].lon, s.elev = jobs[i].elev;
              s.track = *trk[i];
              s.T.assign(g.T, g.T + K), s.d.assign(g.d, g.d + K),
                  s.v.assign(g.v, g.v + K);
              s.f0.assign(g.f0, g.f0 + K), s.score.assign(g.score, g.score + K);
              s.shift = g.shift;
              s.total = g.total;
              valid[i] = 1;
            },
            4096);
        return;
      }
    }
    parallel_for(
        n,
        [&](size_t i, int tid) {
          valid[i] = score(ll[i].first, ll[i].second, *trk[i], *Tg[i], out[i],
                           tb[tid], gb[tid]);
        },
        16);
  }
  // Model curve for overlay plots: (audio time, Hz) for each harmonic.
  std::vector<Curve> curves(size_t k, const Scored& s) const {
    const Event& e0 = ev[0];
    const Event& e = ev[k];
    Vec3 o = geodetic_to_ecef(s.lat, s.lon, s.elev + P.obs_height);
    double dt = s.shift * P.hop;
    int nfr = e.fb1 - e.fb0 + 1;
    std::vector<double> t(nfr);
    std::vector<float> g(nfr);
    for (int i = 0; i < nfr; ++i)
      t[i] = s.T[0] + e.img.time(e.fb0 + i) - e0.tc - dt;
    model_lng(tracks[s.track[k]], o, t.data(), nfr, P.c, g.data());
    std::vector<Curve> out;
    for (double m : e.spec.mult) {
      Curve c{{}, {80, 255, 255}};
      for (int i = 0; i < nfr; ++i)
        if (!std::isnan(g[i]))
          c.tf.push_back(
              {e.img.time(e.fb0 + i), s.f0[k] * m * std::exp((double)g[i])});
      out.push_back(std::move(c));
    }
    return out;
  }
};

std::string flights_json(const std::vector<Flight>& flights,
                         const std::vector<Track>& tracks,
                         const std::vector<uint32_t>& trk) {
  std::string s;
  for (size_t k = 0; k < trk.size(); ++k) {
    if (k) s += " | ";
    s += track_label(flights, tracks[trk[k]]);
  }
  return s;
}

std::string scored_props(const Scored& s, const std::vector<Flight>& flights,
                         const std::vector<Track>& tracks, int rank = -1) {
  char b[256];
  std::snprintf(b, sizeof b, "\"total\":%.4f,\"elev\":%.1f,\"shift_frames\":%d",
                s.total, s.elev, s.shift);
  std::string p = b;
  if (rank >= 0) p += ",\"rank\":" + std::to_string(rank);
  p += ",\"flights\":\"" + json_escape(flights_json(flights, tracks, s.track)) +
       "\"";
  for (size_t k = 0; k < s.track.size(); ++k) {
    std::snprintf(b, sizeof b,
                  ",\"score%zu\":%.4f,\"d%zu\":%.0f,\"v%zu\":%.1f,\"f0_%zu\":%."
                  "2f,\"time%zu\":\"%s\"",
                  k + 1, s.score[k], k + 1, s.d[k], k + 1, s.v[k], k + 1,
                  s.f0[k], k + 1, fmt_time(s.T[k]).c_str());
    p += b;
  }
  return p;
}

}  // namespace

// ----------------------------------------------------------- driver ---
int run_locate(const Audio& audio, const std::vector<Flight>& flights,
               const std::vector<Track>& tracks, const Dem& dem,
               std::vector<Event>& events, const SearchParams& P) {
  (void)audio;
  Timer total;
  fs::create_directories(P.outdir);
  const size_t K = events.size();
  GeoGrid G(P.coarse_step), Gf(P.fine_step);

  // Which flights may have produced the Doppler events: everything that is not
  // a known jet.
  std::vector<char> elig(tracks.size());
  double dmin = 1e300, dmax = -1e300;
  size_t ne = 0;
  for (size_t i = 0; i < tracks.size(); ++i) {
    elig[i] = P.include_jets || tracks[i].cls != ACClass::Jet;
    ne += elig[i];
    dmin = std::min(dmin, tracks[i].t0);
    dmax = std::max(dmax, tracks[i].t1);
  }
  log("tracks: %zu legs, %zu eligible (non-jet), data %s .. %s", tracks.size(),
      ne, fmt_time(dmin).c_str(), fmt_time(dmax).c_str());

  const double a0 = events[0].tc;
  double W0lo = std::max(P.tmin, dmin), W0hi = std::min(P.tmax, dmax);
  if (std::isfinite(P.rec_start)) {
    W0lo =
        std::max(W0lo, P.rec_start + a0 - P.rec_start_tol - events[0].tc_tol);
    W0hi =
        std::min(W0hi, P.rec_start + a0 + P.rec_start_tol + events[0].tc_tol);
  }
  log("event 1 delayed-CPA window: %s .. %s", fmt_time(W0lo).c_str(),
      fmt_time(W0hi).c_str());
  std::vector<double> delta(K), tol(K);
  for (size_t k = 0; k < K; ++k) {
    delta[k] = events[k].tc - a0;
    tol[k] = k ? events[k].tc_tol + events[0].tc_tol + P.offset_slack : 0;
  }

  // ---- A-C, tile by tile under a memory budget
  // -------------------------------------------------- Hits (grid node x flight
  // pass) are by far the largest intermediate: a continental-US day produces
  // hundreds of millions. They only matter until the join, so each spatial tile
  // is generated, joined and silence-filtered on its own and only the (few)
  // surviving candidates are kept. A tile whose hits would exceed the budget is
  // split in four and retried.
  const int nt = num_threads();
  if (gpu::enabled()) {
    Timer tu;
    bool ok = gpu::upload_tracks(tracks) && gpu::upload_dem(dem);
    if (ok) {
      int S = (int)std::ceil(P.dt_sweep / P.hop);
      std::vector<gpu::EventImage> gi;
      for (auto& x : events) {
        gpu::EventImage g{};
        g.himg = x.H.himg.data(), g.hblur3 = x.H.hblur3.data();
        g.nb = x.H.nb, g.E = x.H.E, g.nrH = x.H.nrH, g.nq = x.H.nq;
        for (int i = 0; i < 8; ++i)
          g.thr[i] = x.H.thr[i];
        g.nf = x.img.nf, g.nr = x.img.nr, g.fb0 = x.fb0, g.fb1 = x.fb1;
        g.rlo = (int)std::floor(x.img.row_of(x.f0min)),
        g.rhi = (int)std::ceil(x.img.row_of(x.f0max));
        g.t0 = x.img.t0, g.hop = x.img.hop, g.lnf0 = x.img.lnf0,
        g.dln = x.img.dln, g.tc = x.tc;
        for (double k : x.spec.mult)
          g.mrow.push_back((float)(std::log(k) / x.img.dln));
        gi.push_back(std::move(g));
      }
      ok = gpu::upload_events(gi, S, P.hop, P.c, P.obs_height);
    }
    log("gpu: %s (%.1fs upload)",
        ok ? gpu::info().c_str() : "disabled after upload failure - CPU path",
        tu.sec());
  }
  double Rmax = 0;
  for (auto& e : events)
    Rmax = std::max(Rmax, e.dmax * 1.03);
  const size_t rss0 = current_rss_bytes();
  const size_t limit = P.mem_limit ? P.mem_limit : resolve_mem_limit(0);
  size_t budget = limit > rss0 ? (size_t)((limit - rss0) * 0.6) : 0;
  budget = std::max<size_t>(budget, (size_t)256 << 20);
  const double nodes_disk =
      kPi * Rmax * Rmax / (P.coarse_step * P.coarse_step) +
      4 * Rmax / P.coarse_step + 8;
  const int chunk = (int)std::clamp(
      budget * 0.15 / (nt * nodes_disk * 16 * 2.0), 32.0, 4096.0);
  // hits of all K events coexist during the join; x2 for per-thread parts +
  // merged copy
  const size_t hit_cap = std::max<size_t>(
      100000, (size_t)(budget * 0.8 / (K * sizeof(Hit) * 2.0)));
  log("memory: limit %.1f GB, in use %.1f GB -> %zu hits/event/tile, %d-sample "
      "chunks",
      gb(limit), gb(rss0), hit_cap, chunk);

  std::vector<std::pair<double, double>>
      rel;  // quiet intervals relative to event-1 audio CPA time
  {
    std::vector<double> tcs;
    for (auto& e : events)
      tcs.push_back(e.tc);
    std::sort(tcs.begin(), tcs.end());
    for (size_t i = 0; i + 1 < tcs.size(); ++i)
      if (tcs[i + 1] - tcs[i] > 2 * P.silence_guard)
        rel.push_back(
            {tcs[i] + P.silence_guard - a0, tcs[i + 1] - P.silence_guard - a0});
    for (auto& q : P.quiet)
      rel.push_back({q.first - a0, q.second - a0});
  }
  const bool do_silence = P.silence && !rel.empty();
  const double Rs = std::max(P.silence_radius, P.silence_radius_jet);

  // initial tiles over the extent of eligible traffic
  double ela0 = 90, ela1 = -90, elo0 = 180, elo1 = -180;
  for (size_t i = 0; i < tracks.size(); ++i)
    if (elig[i])
      ela0 = std::min(ela0, tracks[i].minlat),
      ela1 = std::max(ela1, tracks[i].maxlat),
      elo0 = std::min(elo0, tracks[i].minlon),
      elo1 = std::max(elo1, tracks[i].maxlon);
  std::vector<Tile> stack;
  const double td0 = P.tile_deg;
  for (double la = std::floor(ela0 / td0) * td0; la < ela1; la += td0)
    for (double lo = std::floor(elo0 / td0) * td0; lo < elo1; lo += td0)
      stack.push_back({la, la + td0, lo, lo + td0});
  std::reverse(stack.begin(), stack.end());

  // ---- D state: candidates are Doppler-scored as soon as their tile is done
  // (streaming), and only the best location per flight combination is retained
  // -> memory does not grow with candidate count.
  Scorer scorer(tracks, events, dem, P);
  std::vector<std::vector<double>> tbuf(nt);
  std::vector<std::vector<float>> gbuf(nt);
  std::unordered_map<uint64_t, Scored> best;
  GeoJsonWriter dumpD;
  dumpD.open(P.outdir + "/stageD_scored.geojson");
  size_t nD_written = 0, nvalid = 0, njobs = 0, nC = 0;
  const size_t kDumpDCap = 500000;
  double tD = 0;
  auto combo_hash = [](const std::vector<uint32_t>& tr) {
    uint64_t h = 0;
    for (auto t : tr)
      h = mix64(h ^ (t + 0x9e3779b97f4a7c15ULL));
    return h;
  };
  auto score_tile = [&](const std::vector<Candidate>& tc) {
    Timer t;
    struct Job {
      uint64_t fkey;
      uint32_t cand;
    };
    std::vector<Job> jobs;
    std::unordered_set<uint64_t> seen;
    for (uint32_t ci = 0; ci < tc.size(); ++ci) {
      double la, lo;
      G.center(tc[ci].key, la, lo);
      uint64_t th = combo_hash(tc[ci].track);
      Gf.for_each_within(la, lo, P.coarse_step * 0.75,
                         [&](uint64_t fk, double, double) {
                           if (seen.insert(mix64(fk) ^ th).second)
                             jobs.push_back({fk, ci});
                         });
    }
    njobs += jobs.size();
    // Two-deep pipeline: batch b+1 is scored (GPU kernel, or CPU threads) while
    // this thread folds batch b into the best-per-combination table and the
    // stage-D dump.
    const size_t kBatch = 65536;
    struct Buf {
      std::vector<Scored> scored;
      std::vector<char> valid;
      std::vector<std::pair<double, double>> ll;
      std::vector<const std::vector<uint32_t>*> trk;
      std::vector<const std::vector<double>*> Tg;
    };
    Buf buf[2];
    auto launch = [&](size_t b0, Buf& B) {
      size_t n = std::min(kBatch, jobs.size() - b0);
      B.ll.resize(n), B.trk.resize(n), B.Tg.resize(n);
      for (size_t i = 0; i < n; ++i) {
        Gf.center(jobs[b0 + i].fkey, B.ll[i].first, B.ll[i].second);
        B.trk[i] = &tc[jobs[b0 + i].cand].track;
        B.Tg[i] = &tc[jobs[b0 + i].cand].T;
      }
      return std::async(std::launch::async, [&] {
        scorer.score_many(B.ll, B.trk, B.Tg, B.scored, B.valid, tbuf, gbuf);
      });
    };
    std::future<void> fut;
    if (!jobs.empty()) fut = launch(0, buf[0]);
    for (size_t b0 = 0, bi = 0; b0 < jobs.size(); b0 += kBatch, bi ^= 1) {
      fut.get();
      if (b0 + kBatch < jobs.size()) fut = launch(b0 + kBatch, buf[bi ^ 1]);
      Buf& B = buf[bi];
      size_t n = B.ll.size();
      for (size_t i = 0; i < n; ++i) {
        if (!B.valid[i]) continue;
        ++nvalid;
        Scored& sc = B.scored[i];
        if (nD_written < kDumpDCap)
          dumpD.add(sc.lat, sc.lon, scored_props(sc, flights, tracks)),
              ++nD_written;
        uint64_t h = combo_hash(sc.track);
        auto it = best.find(h);
        if (it == best.end())
          best.emplace(h, std::move(sc));
        else if (sc.total > it->second.total)
          it->second = std::move(sc);
      }
    }
    // keep the retained table bounded: prune to the best max_candidates
    // combinations
    if (best.size() > 2 * P.max_candidates) {
      std::vector<std::pair<double, uint64_t>> v;
      for (auto& kv : best)
        v.push_back({kv.second.total, kv.first});
      std::nth_element(v.begin(), v.begin() + P.max_candidates, v.end(),
                       std::greater<>());
      for (size_t i = P.max_candidates; i < v.size(); ++i)
        best.erase(v[i].second);
    }
    tD += t.sec();
  };
  size_t nA_hits = 0, nA_nodes = 0, nB = 0, nB_nodes = 0, nC_nodes = 0,
         ntiles = 0, nsplit = 0, peak_hits = 0;
  const size_t kDumpCap = 500000;
  std::vector<PointFeature> dumpA, dumpB, dumpC;
  bool dumpA_on = true;
  double tA = 0, tB = 0, tC = 0;
  Timer tall;
  std::vector<std::vector<Hit>> hits(K);
  std::vector<uint32_t> sub, sub_sil;
  while (!stack.empty()) {
    Tile T = stack.back();
    stack.pop_back();
    // flights that can matter for this tile
    double mla = Rmax / 110000.0 + 0.03,
           mlo = mla / std::cos(std::min(85.0, std::max(std::fabs(T.la0),
                                                        std::fabs(T.la1))) *
                                kDeg);
    sub.clear();
    for (uint32_t i = 0; i < tracks.size(); ++i) {
      const Track& t = tracks[i];
      if (!elig[i] || t.maxlat < T.la0 - mla || t.minlat > T.la1 + mla ||
          t.maxlon < T.lo0 - mlo || t.minlon > T.lo1 + mlo)
        continue;
      sub.push_back(i);
    }
    if (sub.empty()) continue;
    Timer ta;
    bool ok = true;
    for (size_t k = 0; k < K && ok; ++k)
      ok = gen_hits(events[k], tracks, sub, G, dem, P, W0lo + delta[k] - tol[k],
                    W0hi + delta[k] + tol[k], T, hit_cap, chunk, hits[k]);
    tA += ta.sec();
    if (!ok) {
      for (auto& h : hits)
        std::vector<Hit>().swap(h);
      double hl = (T.la1 - T.la0) / 2, ho = (T.lo1 - T.lo0) / 2;
      if (hl < 0.01)
        die("a %.3f deg tile still exceeds the hit budget (%zu). Raise "
            "--mem-limit or tighten the event ranges.",
            T.la1 - T.la0, hit_cap);
      for (int q = 0; q < 4; ++q) {
        double la = T.la0 + (q / 2) * hl, lo = T.lo0 + (q % 2) * ho;
        stack.push_back({la, la + hl, lo, lo + ho});
      }
      ++nsplit;
      vlog("tile %.2f,%.2f (%.3f deg): over budget, split", T.la0, T.lo0,
           T.la1 - T.la0);
      continue;
    }
    ++ntiles;
    size_t th = 0;
    for (auto& h : hits)
      th += h.size();
    peak_hits = std::max(peak_hits, th);
    nA_hits += hits[0].size();
    {
      std::unordered_set<uint64_t> u;
      if (std::is_sorted(
              hits[0].begin(), hits[0].end(),
              [](const Hit& a, const Hit& b) { return a.key < b.key; })) {
        for (size_t q = 0; q < hits[0].size(); ++q)
          nA_nodes += !q || hits[0][q].key != hits[0][q - 1].key;
        if (dumpA_on && nA_nodes <= 300000)
          for (auto& h : hits[0])
            u.insert(h.key);
        else
          dumpA_on = false;
      } else {
        for (auto& h : hits[0])
          u.insert(h.key);
        nA_nodes += u.size();
      }
      if (dumpA_on) {
        for (auto& h : hits[0])
          if (u.erase(h.key)) {
            double la, lo;
            G.center(h.key, la, lo);
            dumpA.push_back({la, lo,
                             "\"d\":" + std::to_string((int)h.d) +
                                 ",\"time\":\"" + fmt_time(h.T) + "\""});
          }
        if (dumpA.size() > 300000)
          dumpA_on = false, std::vector<PointFeature>().swap(dumpA);
      }
    }
    // ---- B: time-offset join (within the tile)
    Timer tb;
    auto by_key = [](const Hit& a, const Hit& b) { return a.key < b.key; };
    for (size_t k = 1; k < K; ++k)
      if (!std::is_sorted(hits[k].begin(), hits[k].end(), by_key))
        std::sort(hits[k].begin(), hits[k].end(), by_key);
    // parallel over contiguous chunks of hits[0]; per-chunk outputs are
    // concatenated in order (deterministic)
    std::vector<Candidate> tc;
    {
      const size_t nh0 = hits[0].size(),
                   nchunk = std::min<size_t>(nh0, (size_t)nt * 8);
      std::vector<std::vector<Candidate>> part(nchunk);
      parallel_for(
          nchunk,
          [&](size_t ci, int) {
            std::vector<std::array<const Hit*, 8>> opts(K);
            std::vector<int> nopt(K), idx(K);
            auto& out = part[ci];
            for (size_t q = nh0 * ci / nchunk, qe = nh0 * (ci + 1) / nchunk;
                 q < qe; ++q) {
              const Hit& h0 = hits[0][q];
              bool okj = true;
              for (size_t k = 1; k < K && okj; ++k) {
                nopt[k] = 0;
                auto it = std::lower_bound(
                    hits[k].begin(), hits[k].end(), h0.key,
                    [](const Hit& a, uint64_t key) { return a.key < key; });
                for (; it != hits[k].end() && it->key == h0.key; ++it) {
                  if (std::fabs(it->T - h0.T - delta[k]) > tol[k]) continue;
                  if (it->track == h0.track && std::fabs(it->T - h0.T) < 30)
                    continue;
                  if (nopt[k] < 8) opts[k][nopt[k]++] = &*it;
                }
                okj = nopt[k] > 0;
              }
              if (!okj) continue;
              std::fill(idx.begin(), idx.end(), 0);
              for (;;) {
                Candidate c{h0.key, {h0.track}, {h0.T}};
                c.track.reserve(K), c.T.reserve(K);
                for (size_t k = 1; k < K; ++k) {
                  c.track.push_back(opts[k][idx[k]]->track);
                  c.T.push_back(opts[k][idx[k]]->T);
                }
                out.push_back(std::move(c));
                size_t k = 1;
                while (k < K && ++idx[k] == nopt[k])
                  idx[k++] = 0;
                if (k >= K) break;
              }
            }
          },
          1);
      size_t tot = 0;
      for (auto& v : part)
        tot += v.size();
      tc.reserve(tot);
      for (auto& v : part)
        for (auto& c : v)
          tc.push_back(std::move(c));
    }
    // Keep the hit buffers' pages for the next tile (re-faulting ~0.5 GB per
    // tile costs ~0.2 s) unless they are a large share of the memory limit;
    // then free before the silence index is built.
    for (auto& h : hits) {
      if (h.capacity() * sizeof(Hit) > limit / 8)
        std::vector<Hit>().swap(h);
      else
        h.clear();
    }
    auto add_dump = [&](std::vector<PointFeature>& d,
                        const std::vector<Candidate>& v) -> size_t {
      std::unordered_set<uint64_t> seen;
      for (auto& c : v)
        if (seen.insert(c.key).second && d.size() < kDumpCap) {
          double la, lo;
          G.center(c.key, la, lo);
          d.push_back({la, lo,
                       "\"flights\":\"" +
                           json_escape(flights_json(flights, tracks, c.track)) +
                           "\""});
        }
      return seen.size();
    };
    nB += tc.size();
    nB_nodes += add_dump(dumpB, tc);
    tB += tb.sec();
    // ---- C: silence filter (index only over flights near this tile, minutes
    // these candidates need)
    if (do_silence && !tc.empty()) {
      Timer tcs;
      MinuteSet minutes;
      {
        double tlo = 1e300, thi = -1e300, rlo_ = 1e300, rhi_ = -1e300;
        for (auto& c : tc)
          tlo = std::min(tlo, c.T[0]), thi = std::max(thi, c.T[0]);
        for (auto& r : rel)
          rlo_ = std::min(rlo_, r.first), rhi_ = std::max(rhi_, r.second);
        minutes.build((int64_t)std::floor((tlo + rlo_ - 30) / 60) - 1,
                      (int64_t)std::floor((thi + rhi_ + 60) / 60) + 1);
        for (auto& c : tc)
          for (auto& r : rel)
            for (double x = c.T[0] + r.first - 30; x <= c.T[0] + r.second + 60;
                 x += 60)
              minutes.insert((int64_t)std::floor(x / 60));
      }
      double sla = Rs / 110000.0 + 0.05,
             slo = sla / std::cos(std::min(85.0, std::fabs(T.la1)) * kDeg);
      sub_sil.clear();
      for (uint32_t i = 0; i < tracks.size(); ++i) {
        const Track& t = tracks[i];
        if (t.maxlat < T.la0 - sla || t.minlat > T.la1 + sla ||
            t.maxlon < T.lo0 - slo || t.minlon > T.lo1 + slo)
          continue;
        sub_sil.push_back(i);
      }
      double tC_pre = tcs.sec();
      SilenceIndex idx(tracks, sub_sil, minutes, 4.0, T.la0, T.la1, T.lo0,
                       T.lo1);
      double tC_idx = tcs.sec();
      std::vector<char> keep(tc.size(), 0);
      parallel_for(
          tc.size(),
          [&](size_t i, int) {
            const Candidate& c = tc[i];
            double la, lo;
            G.center(c.key, la, lo);
            Vec3 o = geodetic_to_ecef(la, lo, dem.elev(la, lo) + P.obs_height);
            auto excl = [&](uint32_t tr, double trx) {
              for (size_t k = 0; k < K; ++k)
                if (c.track[k] == tr &&
                    std::fabs(trx - c.T[k]) < P.silence_guard + 90)
                  return true;
              return false;
            };
            for (auto& r : rel)
              if (idx.heard(o, la, lo, c.T[0] + r.first, c.T[0] + r.second,
                            P.silence_radius, P.silence_radius_jet, P.c, excl))
                return;
            keep[i] = 1;
          },
          64);
      std::vector<Candidate> k2;
      for (size_t i = 0; i < tc.size(); ++i)
        if (keep[i]) k2.push_back(std::move(tc[i]));
      tc.swap(k2);
      vlog("silence filter: minutes+subset %.2fs  index %.2fs  queries %.2fs",
           tC_pre, tC_idx - tC_pre, tcs.sec() - tC_idx);
      tC += tcs.sec();
    }
    nC_nodes += add_dump(dumpC, tc);
    nC += tc.size();
    if (!tc.empty()) score_tile(tc);
    if (g_verbose)
      std::fprintf(stderr,
                   "  ... tiles %zu done (%zu split), %zu queued | %zu "
                   "candidates, %zu scored | rss %.1f GB   \r",
                   ntiles, nsplit, stack.size(), nC, nvalid,
                   gb(current_rss_bytes()));
  }
  if (g_verbose) std::fprintf(stderr, "\n");
  log("[A] event-1 CPA hits: %zu at %zu grid nodes over %zu tiles (%zu splits, "
      "peak %zu hits in one tile)  %.1fs",
      nA_hits, nA_nodes, ntiles, nsplit, peak_hits, tA);
  log("[B] time-offset filter: %zu candidates at %zu grid nodes  %.1fs", nB,
      nB_nodes, tB);
  if (do_silence)
    log("[C] silence filter: %zu candidates at %zu grid nodes remain  %.1fs",
        nC, nC_nodes, tC);
  else
    log("[C] silence filter: skipped (%s)",
        P.silence ? "no quiet interval" : "disabled");
  log("[D] Doppler modelling: %zu (node, flight-combo) jobs on %.0f m grid, "
      "%zu scored  %.1fs",
      njobs, P.fine_step, nvalid, tD);
  log("stages A-D: %.1fs", tall.sec());
  dumpD.close();
  if (nD_written >= kDumpDCap)
    log("note: stageD_scored.geojson truncated to the first %zu nodes",
        kDumpDCap);
  if (!dumpA.empty())
    write_points_geojson(P.outdir + "/stageA_event1_hits.geojson", dumpA);
  write_points_geojson(P.outdir + "/stageB_offset.geojson", dumpB);
  if (do_silence)
    write_points_geojson(P.outdir + "/stageC_silence.geojson", dumpC);
  std::vector<PointFeature>().swap(dumpA),
      std::vector<PointFeature>().swap(dumpB),
      std::vector<PointFeature>().swap(dumpC);
  if (!nC) {
    log("no candidates survived. Loosen ranges "
        "(--v-margin/--d-margin/--offset-slack), check times, or add data.");
    return 2;
  }
  if (!nvalid) {
    log("no location could be modelled");
    return 2;
  }

  // ---- E: best per flight combination
  std::vector<Scored> ranked;
  ranked.reserve(best.size());
  for (auto& kv : best)
    ranked.push_back(std::move(kv.second));
  std::unordered_map<uint64_t, Scored>().swap(best);
  std::sort(ranked.begin(), ranked.end(),
            [](const Scored& a, const Scored& b) { return a.total > b.total; });
  size_t ntop = std::min<size_t>(P.top, ranked.size());

  std::printf(
      "\n=== ranked locations (best per flight combination; higher is better) "
      "===\n");
  std::printf("%-4s %-11s %-12s %-8s", "rank", "lat", "lon", "total");
  for (size_t k = 0; k < K; ++k)
    std::printf(" ev%zu   ", k + 1);
  std::printf(" flights\n");
  for (size_t i = 0; i < ntop; ++i) {
    const Scored& s = ranked[i];
    std::printf("%-4zu %-11.6f %-12.6f %-8.3f", i + 1, s.lat, s.lon, s.total);
    for (size_t k = 0; k < K; ++k)
      std::printf(" %-6.3f", s.score[k]);
    std::printf(" %s\n", flights_json(flights, tracks, s.track).c_str());
  }
  if (ranked.size() > 1)
    std::printf("margin #1 over #2: %.3f (%.1f%%)\n",
                ranked[0].total - ranked[1].total,
                100.0 * (ranked[0].total - ranked[1].total) /
                    std::max(1e-9, std::fabs(ranked[1].total)));
  {
    std::vector<PointFeature> pf;
    for (size_t i = 0; i < ranked.size(); ++i)
      pf.push_back({ranked[i].lat, ranked[i].lon,
                    scored_props(ranked[i], flights, tracks, (int)i + 1)});
    write_points_geojson(P.outdir + "/ranked.geojson", pf);
    pf.resize(std::min<size_t>(1, pf.size()));
    write_points_geojson(P.outdir + "/best.geojson", pf);
  }

  // Overlays: spectrogram + modelled curve for the top results (the "does it
  // cover the curve" check).
  for (size_t i = 0; i < std::min<size_t>(P.overlays, ranked.size()); ++i) {
    for (size_t k = 0; k < K; ++k) {
      const Event& e = events[k];
      Image img = render_spectrogram(e.white, e.spec.t0 - 5, e.spec.t1 + 5,
                                     e.spec.fmin * 0.85, e.spec.fmax * 1.15,
                                     scorer.curves(k, ranked[i]), 2);
      write_png(P.outdir + "/overlay_rank" + std::to_string(i + 1) + "_event" +
                    std::to_string(k + 1) + ".png",
                img);
    }
  }

  // Refinement heatmap around the winners using only their flights (dense
  // grid).
  std::vector<std::string> heat_files;
  std::vector<std::vector<double>> heat_vals;
  for (size_t i = 0; i < std::min<size_t>(P.refine_top, ranked.size()); ++i) {
    Timer tr;
    GeoGrid Gr(P.refine_step);
    std::vector<uint64_t> keys;
    Gr.for_each_within(ranked[i].lat, ranked[i].lon, P.refine_radius,
                       [&](uint64_t k, double, double) { keys.push_back(k); });
    std::vector<Scored> rs;
    std::vector<char> ok;
    {
      std::vector<std::pair<double, double>> ll(keys.size());
      for (size_t j = 0; j < keys.size(); ++j)
        Gr.center(keys[j], ll[j].first, ll[j].second);
      std::vector<const std::vector<uint32_t>*> trk(keys.size(),
                                                    &ranked[i].track);
      std::vector<const std::vector<double>*> Tg(keys.size(), &ranked[i].T);
      scorer.score_many(ll, trk, Tg, rs, ok, tbuf, gbuf);
    }
    std::vector<PointFeature> pf;
    heat_vals.emplace_back();
    size_t bi = 0;
    double bv = -1e30;
    for (size_t j = 0; j < keys.size(); ++j)
      if (ok[j]) {
        heat_vals.back().push_back(rs[j].total);
        pf.push_back(
            {rs[j].lat, rs[j].lon, scored_props(rs[j], flights, tracks)});
        if (rs[j].total > bv) bv = rs[j].total, bi = j;
      }
    std::string f = "refine_rank" + std::to_string(i + 1) + ".geojson";
    write_points_geojson(P.outdir + "/" + f, pf);
    heat_files.push_back(f);
    if (bv > ranked[i].total) {
      log("[E] refine rank %zu: %zu nodes @%.0f m, peak %.4f at %.6f, %.6f "
          "(%.1fs)",
          i + 1, pf.size(), P.refine_step, bv, rs[bi].lat, rs[bi].lon,
          tr.sec());
      if (i == 0) ranked[0] = rs[bi];
    } else {
      log("[E] refine rank %zu: %zu nodes @%.0f m (%.1fs)", i + 1, pf.size(),
          P.refine_step, tr.sec());
    }
  }

  // Flight paths of the top combinations (clipped to +-10 min around the
  // events).
  {
    std::vector<LineFeature> lf;
    std::unordered_set<uint32_t> seen;
    for (size_t i = 0; i < std::min<size_t>(5, ranked.size()); ++i)
      for (size_t k = 0; k < K; ++k) {
        uint32_t ti = ranked[i].track[k];
        if (!seen.insert(ti).second) continue;
        const Track& tr = tracks[ti];
        LineFeature L;
        int hint = -1;
        Vec3 p, v;
        for (double t = std::max(tr.t0, ranked[i].T[k] - 600);
             t <= std::min(tr.t1, ranked[i].T[k] + 600); t += 2) {
          tr.eval(t, p, v, &hint);
          double la, lo, h;
          ecef_to_geodetic(p, la, lo, h);
          L.lonlat.push_back({lo, la});
        }
        if (L.lonlat.size() < 2) continue;
        L.props = "\"flight\":\"" + json_escape(track_label(flights, tr)) +
                  "\",\"event\":" + std::to_string(k + 1) +
                  ",\"rank\":" + std::to_string(i + 1) + ",\"class\":\"" +
                  class_name(tr.cls) + "\"";
        lf.push_back(std::move(L));
      }
    write_lines_geojson(P.outdir + "/flights.geojson", lf);
  }

  const Scored& B = ranked[0];
  std::printf("\nBEST: %.6f, %.6f  (elev %.0f m)  total score %.3f\n", B.lat,
              B.lon, B.elev, B.total);
  for (size_t k = 0; k < K; ++k)
    std::printf(
        "  event %zu: %s  CPA %.0f m @ %.1f m/s  heard %s  f0 %.2f Hz  score "
        "%.3f\n",
        k + 1, track_label(flights, tracks[B.track[k]]).c_str(), B.d[k], B.v[k],
        fmt_time(B.T[k]).c_str(), B.f0[k], B.score[k]);
  double rec0 = B.T[0] - events[0].tc - B.shift * P.hop;
  std::printf("  implied recording start (audio t=0): %s\n",
              fmt_time(rec0).c_str());
  std::printf(
      "  "
      "https://www.google.com/maps/"
      "@?api=1&map_action=pano&viewpoint=%.6f,%.6f\n",
      B.lat, B.lon);
  std::printf("  https://www.google.com/maps/search/?api=1&query=%.6f,%.6f\n",
              B.lat, B.lon);

  // result.json
  {
    std::ofstream f(P.outdir + "/result.json");
    f << "{\n  \"best\": {\"lat\": " << B.lat << ", \"lon\": " << B.lon
      << ", \"total\": " << B.total << ", \"recording_start\": \""
      << fmt_time(rec0) << "\"},\n  \"events\": [\n";
    for (size_t k = 0; k < K; ++k) {
      const Event& e = events[k];
      char b[1024];
      std::snprintf(
          b, sizeof b,
          "    {\"box\": [%.2f, %.2f, %.1f, %.1f], \"fit\": {\"v\": %.2f, "
          "\"d\": %.1f, \"tc\": %.3f, "
          "\"f0\": %.3f}, \"ranges\": {\"v\": [%.2f, %.2f], \"d\": [%.1f, "
          "%.1f], \"tc_tol\": %.2f, "
          "\"f0\": [%.3f, %.3f]}, \"flight\": \"%s\", \"heard\": \"%s\", "
          "\"score\": %.4f}%s\n",
          e.spec.t0, e.spec.t1, e.spec.fmin, e.spec.fmax, e.fit.v, e.fit.d,
          e.fit.tc, e.fit.f0, e.vmin, e.vmax, e.dmin, e.dmax, e.tc_tol, e.f0min,
          e.f0max,
          json_escape(track_label(flights, tracks[B.track[k]])).c_str(),
          fmt_time(B.T[k]).c_str(), B.score[k], k + 1 < K ? "," : "");
      f << b;
    }
    f << "  ],\n  \"ranked\": [\n";
    for (size_t i = 0; i < ntop; ++i) {
      char b[512];
      std::snprintf(
          b, sizeof b,
          "    {\"rank\": %zu, \"lat\": %.6f, \"lon\": %.6f, \"total\": %.4f, "
          "\"flights\": \"%s\"}%s\n",
          i + 1, ranked[i].lat, ranked[i].lon, ranked[i].total,
          json_escape(flights_json(flights, tracks, ranked[i].track)).c_str(),
          i + 1 < ntop ? "," : "");
      f << b;
    }
    f << "  ]\n}\n";
  }

  // QGIS project
  ProjectSpec proj;
  proj.title = "audio2loc result";
  double span = std::max(0.02, P.refine_radius / 111000.0 * 2.5);
  proj.minlat = B.lat - span;
  proj.maxlat = B.lat + span;
  proj.minlon = B.lon - span / std::cos(B.lat * kDeg);
  proj.maxlon = B.lon + span / std::cos(B.lat * kDeg);
  {
    LayerSpec L;
    L.name = "best location";
    L.file = "best.geojson";
    L.color = "255,0,0";
    L.size = 6;
    L.marker = "star";
    proj.layers.push_back(L);
  }
  {
    LayerSpec L;
    L.name = "ranked (best per flight combo)";
    L.file = "ranked.geojson";
    L.style = LayerSpec::Graduated;
    L.field = "total";
    std::vector<double> v;
    for (auto& r : ranked)
      v.push_back(r.total);
    L.breaks = quantile_breaks(v, 5);
    L.size = 3;
    proj.layers.push_back(L);
  }
  for (size_t i = 0; i < heat_files.size(); ++i) {
    LayerSpec L;
    L.name = "refine heatmap rank " + std::to_string(i + 1);
    L.file = heat_files[i];
    L.style = LayerSpec::Graduated;
    L.field = "total";
    L.breaks = quantile_breaks(heat_vals[i], 7);
    L.size = 1.2;
    L.visible = i == 0;
    proj.layers.push_back(L);
  }
  {
    LayerSpec L;
    L.name = "flight paths (top 5)";
    L.file = "flights.geojson";
    L.kind = LayerSpec::Lines;
    L.color = "255,0,200";
    L.size = 0.6;
    proj.layers.push_back(L);
  }
  {
    LayerSpec L;
    L.name = "D: Doppler-scored nodes";
    L.file = "stageD_scored.geojson";
    L.style = LayerSpec::Heatmap;
    L.field = "total";
    L.visible = false;
    proj.layers.push_back(L);
  }
  for (auto [n, f, col] :
       std::vector<std::tuple<const char*, const char*, const char*>>{
           {"C: after silence filter", "stageC_silence.geojson", "40,160,255"},
           {"B: after time-offset filter", "stageB_offset.geojson",
            "120,120,120"},
           {"A: event-1 CPA hits", "stageA_event1_hits.geojson",
            "200,200,200"}}) {
    if (!fs::exists(P.outdir + "/" + f)) continue;
    LayerSpec L;
    L.name = n;
    L.file = f;
    L.color = col;
    L.size = 1.0;
    L.visible = false;
    proj.layers.push_back(L);
  }
  write_qgs_project(P.outdir + "/audio2loc.qgs", proj);
  log("wrote QGIS project %s/audio2loc.qgs", P.outdir.c_str());
  if (!P.qgis_native_out.empty()) {
    if (write_qgis_native(P.qgis_native_out, P.outdir, proj))
      log("wrote %s via native QGIS API", P.qgis_native_out.c_str());
    else
      log("native QGIS export unavailable in this build (configure with "
          "-DAUDIO2LOC_WITH_QGIS=ON)");
  }
  log("done in %.1fs -> %s/  (peak memory %.2f GB of %.1f GB limit)",
      total.sec(), P.outdir.c_str(), gb(peak_rss_bytes()), gb(P.mem_limit));
  return 0;
}

}  // namespace pl
