// Optional CUDA acceleration (build with -DAUDIO2LOC_CUDA=ON). Every entry
// point returns false when the GPU path is unavailable so callers fall back to
// the CPU implementation, which stays the reference: GPU kernels mirror the CPU
// algorithms step for step.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "adsb.h"
#include "dem.h"

namespace pl::gpu {

constexpr int kMaxEvents = 6;

bool compiled();                // built with CUDA support
bool init(bool want, int dev);  // select device; returns enabled()
bool enabled();
std::string info();

// ------------------------------------------------------------------ DSP ---
// Hann-windowed STFT magnitude in dB, identical layout/normalisation to
// compute_spectrogram().
bool stft_db(const float* x, size_t n, double sr, int wl, int nfft, double t0,
             double hop, int nframes, int bin0, int nbins, float* out);
// Robust straight-line Doppler fit cost over a (v, d, tc) grid (see doppler.cpp
// fit_cost).
bool fit_costs(const std::vector<double>& t, const std::vector<double>& lnf,
               const std::vector<double>& w, const std::vector<double>& vs,
               const std::vector<double>& ds, const std::vector<double>& tcs,
               double c, int side, std::vector<float>& cost,
               std::vector<float>& lnf0);

// --------------------------------------------------------------- search ---
bool upload_tracks(const std::vector<Track>& tracks);
bool upload_dem(const Dem& dem);

struct HitParams {
  double vmin, vmax, dmin, dmax, Tlo, Thi, c, dt, obs_h, step;
  double la0, la1, lo0, lo1;  // owning tile
};
struct SampleRange {  // samples k in [ka, kb] at t = a + k*dt of one track
  uint32_t track;
  int32_t ka, kb;
  double a;
};
struct GHit {  // same layout as search.cpp Hit
  uint64_t key;
  uint32_t track;
  float d, v;
  double T;
};
// 1 = ok, 0 = more than `cap` hits (caller splits the tile), -1 = GPU failure
// (fall back to CPU).
int gen_hits(const HitParams& hp, const std::vector<SampleRange>& ranges,
             size_t cap, std::vector<GHit>& out);

struct EventImage {
  const float*
      himg;  // harmonic-summed tables over the box frames (see Event::Harm)
  const float* hblur3;  // [(j*nb + b)*3 + p][nq]
  int nb, E, nrH, nq;
  float thr[8];
  int nf, nr, fb0, fb1, rlo, rhi;
  double t0, hop, lnf0, dln, tc;
  std::vector<float> mrow;  // harmonic row offsets
};
bool upload_events(const std::vector<EventImage>& ev, int S, double hop,
                   double c, double obs_h);

struct ScoreJob {
  double lat, lon, elev;
  uint32_t track[kMaxEvents];
  double Tguess[kMaxEvents];
};
struct ScoreOut {
  int valid, shift;
  float total;
  float score[kMaxEvents];
  double T[kMaxEvents], d[kMaxEvents], v[kMaxEvents], f0[kMaxEvents];
};
bool score(const std::vector<ScoreJob>& jobs, int K,
           std::vector<ScoreOut>& out);

}  // namespace pl::gpu
