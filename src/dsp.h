// Audio I/O, FFT, spectrograms and the log-frequency scoring images.
#pragma once
#include <complex>
#include <string>
#include <vector>

namespace pl {

struct Audio {
  double sr = 0;
  std::vector<float> x;  // mono
  double duration() const { return sr > 0 ? x.size() / sr : 0; }
};

// Loads WAV natively; anything else (mp4, mkv, mp3, m4a, ...) is decoded through ffmpeg.
// Output is mono, low-pass filtered and decimated to roughly target_sr (>= target_sr).
Audio load_audio(const std::string& path, double target_sr, const std::string& ffmpeg = "ffmpeg");
void write_wav16(const std::string& path, const Audio& a);
Audio decimate(const Audio& a, double target_sr);

class FFT {
 public:
  explicit FFT(int n);
  int size() const { return n_; }
  void forward(std::complex<float>* d) const;

 private:
  int n_;
  std::vector<int> rev_;
  std::vector<std::complex<float>> tw_;
};

// Magnitude spectrogram in dB. Frame f is centred at t0 + f*hop (audio seconds).
struct Spectrogram {
  double t0 = 0, hop = 0, win = 0, df = 0;
  int nframes = 0, nbins = 0, bin0 = 0;  // bins [bin0, bin0+nbins) of the FFT
  std::vector<float> v;                  // nframes * nbins, frame-major
  float at(int f, int b) const { return v[(size_t)f * nbins + b]; }
  float& at(int f, int b) { return v[(size_t)f * nbins + b]; }
  double freq(int b) const { return (bin0 + b) * df; }
  double time(int f) const { return t0 + f * hop; }
};

Spectrogram compute_spectrogram(const Audio& a, double win_s, double hop_s, double t_begin, double t_end,
                                double fmin, double fmax, int zero_pad = 2);
// Per-frequency background removal: z = (dB - median_t) / (1.4826*MAD_t), clipped to [0, clip].
void whiten(Spectrogram& s, float clip = 8.f);

// Spectrogram resampled onto a log-frequency axis: row r <-> f = exp(lnf0 + r*dln).
// In this domain a Doppler curve f(t) = f0 * g(t) shifts vertically by an integer when f0 changes,
// and harmonics m*f0 are a constant row offset ln(m)/dln.
struct LogImage {
  double t0 = 0, hop = 0, lnf0 = 0, dln = 0;
  int nf = 0, nr = 0;
  std::vector<float> v;  // nf * nr
  float at(int f, int r) const {
    if ((unsigned)f >= (unsigned)nf || (unsigned)r >= (unsigned)nr) return 0.f;
    return v[(size_t)f * nr + r];
  }
  int frame_of(double t) const { return (int)std::lround((t - t0) / hop); }
  double time(int f) const { return t0 + f * hop; }
  double row_of(double f) const { return (std::log(f) - lnf0) / dln; }
  double freq(double r) const { return std::exp(lnf0 + r * dln); }
  LogImage blurred(int rf, int rr) const;  // box blur (separable), radius in frames / rows
};

LogImage to_log_image(const Spectrogram& whitened, double fmin, double fmax, double dln);

}  // namespace pl
