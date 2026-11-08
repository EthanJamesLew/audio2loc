// The geolocation search.
#pragma once
#include <string>
#include <utility>
#include <vector>

#include "adsb.h"
#include "dem.h"
#include "doppler.h"
#include "dsp.h"
#include "util.h"

namespace pl {

// What the user tells us about one Doppler curve in the recording.
struct EventSpec {
  double t0 = 0, t1 = 0, fmin = 0, fmax = 0;  // box: audio seconds, Hz (around ONE tone of the curve)
  std::vector<double> mult{1.0};              // harmonics (multiples of the boxed tone) used in scoring
  int side = 0;                               // fit: 0 both sides, -1 approach only, +1 recede only
  bool ridge_score = false;                   // score vs traced ridge instead of raw brightness (faint curves)
  std::string trace;                          // optional CSV "time_s,freq_hz" hand-traced curve (fit + ridge score)
  // Optional manual overrides of the measured ranges (NaN = use measurement).
  double vmin = kNaN, vmax = kNaN, dmin = kNaN, dmax = kNaN, tc = kNaN, tc_tol = kNaN, f0min = kNaN, f0max = kNaN;
};

struct Event {
  EventSpec spec;
  DopplerFit fit;
  Ridge ridge;
  double vmin, vmax, dmin, dmax, tc, tc_tol, f0min, f0max;  // final search ranges
  Spectrogram white;   // whitened linear spectrogram around the box (for overlays)
  LogImage img, blur;  // scoring images
  // Harmonic-summed scoring tables over the box frames [fb0, fb1] (built by build_harmonic_tables).
  // eval_curve sums img(f, lrint(x + mrow_m) + r) over harmonics m. Writing x = fl + phi (fl = floor x),
  // lrint(phi + mrow_m) only changes where phi crosses 0.5 - frac(mrow_m), so phi falls in one of nb <= nm+1
  // intervals and each interval has a fixed integer offset pattern off[b][m]. Pre-summing the harmonics per
  // pattern makes every model point ONE table read instead of nm, exactly.
  struct Harm {
    int nb = 1, E = 0, nrH = 0, nq = 0, nfr = 0;
    float thr[8] = {};  // sorted interval boundaries in [0,1) (nb-1 of them)
    // himg[(j*nb + b)*nrH + rr] = sum_m img(fb0+j, rr - E + off[b][m])
    // hblur3[((j*nb + b)*3 + p)*nq + q] = same over blur, row rr = 3q + p (row phase split for the step-3 grid)
    std::vector<float> himg, hblur3;
    int variant(float x, int& fl) const {
      float f = std::floor(x);
      fl = (int)f;
      float phi = x - f;
      int b = 0;
      while (b < nb - 1 && phi >= thr[b]) ++b;
      return b;
    }
  } H;
  int fb0 = 0, fb1 = 0;  // box frame range inside img
};

struct SearchParams {
  double c = 343.0;
  // spectrogram
  double win = 1.5, hop = 0.1, dln = 0.0015, target_sr = 4000;
  // measurement margins (the video: "take a pretty generous range")
  double v_margin = 0.10, d_margin = 0.30, tc_margin = 2.0, f0_margin = 0.01;
  // grids
  double coarse_step = 200, fine_step = 100, sample_dt = 2, obs_height = 2;
  double offset_slack = 5;  // extra +-s allowed on inter-event offsets in the coarse join
  // silence filter
  bool silence = true;
  double silence_radius = 2000, silence_radius_jet = 4000, silence_guard = 60;
  std::vector<std::pair<double, double>> quiet;  // extra audio intervals with no aircraft audible
  // time constraints
  double tmin = -1e300, tmax = 1e300;            // data / search window (absolute)
  double rec_start = kNaN, rec_start_tol = 300;  // absolute time of audio t=0 (optional)
  bool include_jets = false;
  // fine scoring
  double dt_sweep = 3.0;  // +- seconds of anchor alignment sweep
  int top = 25, overlays = 3, refine_top = 1;
  double refine_step = 25, refine_radius = 1500;
  size_t max_candidates = 400000;
  size_t mem_limit = 0;    // bytes; 0 = 50% of currently available RAM
  double tile_deg = 2.0;   // initial tile size for stages A-C (auto-split under memory pressure)
  std::string outdir = "planeloc_out";
  std::string qgis_native_out;  // if non-empty and built with QGIS: write .qgz via QGIS API
};

// Measure each event from the audio (ridge + straight-line Doppler fit) and build scoring images.
std::vector<Event> prepare_events(const Audio& a, const std::vector<EventSpec>& specs, const SearchParams& P);
void print_events(const std::vector<Event>& ev);

struct Scored {
  double lat = 0, lon = 0, elev = 0;
  std::vector<uint32_t> track;
  std::vector<double> T, d, v, f0, score;  // per event: delayed-CPA abs time, slant dist, speed, f0, score
  int shift = 0;                           // anchor alignment shift (frames)
  double total = -1e30;
};

int run_locate(const Audio& audio, const std::vector<Flight>& flights, const std::vector<Track>& tracks,
               const Dem& dem, std::vector<Event>& events, const SearchParams& P);

}  // namespace pl
