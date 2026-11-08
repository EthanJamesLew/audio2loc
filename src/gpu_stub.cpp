// CPU-only build: every GPU entry point reports "unavailable" so callers use
// the CPU path.
#include "gpu.h"

namespace pl::gpu {
bool compiled() {
  return false;
}
bool init(bool, int) {
  return false;
}
bool enabled() {
  return false;
}
std::string info() {
  return "not compiled (configure with -DAUDIO2LOC_CUDA=ON)";
}
bool stft_db(const float*, size_t, double, int, int, double, double, int, int,
             int, float*) {
  return false;
}
bool fit_costs(const std::vector<double>&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<double>&,
               const std::vector<double>&, const std::vector<double>&, double,
               int, std::vector<float>&, std::vector<float>&) {
  return false;
}
bool upload_tracks(const std::vector<Track>&) {
  return false;
}
bool upload_dem(const Dem&) {
  return false;
}
int gen_hits(const HitParams&, const std::vector<SampleRange>&, size_t,
             std::vector<GHit>&) {
  return -1;
}
bool upload_events(const std::vector<EventImage>&, int, double, double,
                   double) {
  return false;
}
bool score(const std::vector<ScoreJob>&, int, std::vector<ScoreOut>&) {
  return false;
}
}  // namespace pl::gpu
