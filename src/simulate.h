#pragma once
#include <cstdint>
#include <string>

namespace pl {

struct SimOptions {
  std::string prefix = "sim";
  std::string start = "2025-09-08T12:00:00Z";
  double hours = 6;
  double lat = 40.85, lon = -73.45,
         region = 1.2;  // centre and square side (deg)
  int n_prop = 300, n_jet = 150, n_decoy = 30;
  double gap = 920;   // seconds between the two flyovers (as in the video)
  double lead = 240;  // recording margin before event 1 / after event 2
  double sr = 8000, noise = 0.08, amp_ref = 0.08, c = 343.0;
  double box = 55;          // half-width (s) of suggested event boxes
  double max_slant = 1200;  // real-traffic mode: max CPA slant distance of the
                            // two event passes (m)
  uint64_t seed = 7;
};

int run_simulate(const SimOptions& o);
// Scenario built on real ADS-B (a .plb/.csv/trace source): real tracks, hidden
// observer, synthetic audio.
int run_simulate_real(const SimOptions& o, const std::string& adsb);

}  // namespace pl
