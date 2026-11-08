#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <random>
#include <vector>

#include "dsp.h"
#include "geo.h"

using Catch::Approx;
using namespace pl;

namespace {

Audio tone(double sr, double seconds, double hz, double amp = 0.5) {
  Audio a;
  a.sr = sr;
  a.x.resize(static_cast<size_t>(sr * seconds));
  for (size_t i = 0; i < a.x.size(); ++i)
    a.x[i] = static_cast<float>(amp * std::sin(2 * kPi * hz * i / sr));
  return a;
}

}  // namespace

TEST_CASE("FFT matches a direct DFT", "[dsp]") {
  constexpr int kN = 64;
  std::mt19937 rng(1);
  std::uniform_real_distribution<float> u(-1.f, 1.f);
  std::vector<std::complex<float>> x(kN), y;
  for (auto& v : x)
    v = {u(rng), u(rng)};
  y = x;
  FFT(kN).forward(y.data());
  for (int k = 0; k < kN; ++k) {
    std::complex<double> s = 0;
    for (int n = 0; n < kN; ++n)
      s += std::complex<double>(x[n]) * std::polar(1.0, -2 * kPi * k * n / kN);
    CHECK(y[k].real() == Approx(s.real()).margin(1e-4));
    CHECK(y[k].imag() == Approx(s.imag()).margin(1e-4));
  }
}

TEST_CASE("a spectrogram peaks at the tone's frequency", "[dsp]") {
  const Audio a = tone(8000, 4.0, 137.0);
  const Spectrogram s = compute_spectrogram(a, 1.0, 0.25, 0, 4, 40, 400);
  REQUIRE(s.nframes > 4);
  const int f = s.nframes / 2;
  int best = 0;
  for (int b = 1; b < s.nbins; ++b)
    if (s.at(f, b) > s.at(f, best)) best = b;
  CHECK(s.freq(best) == Approx(137.0).margin(s.df));
}

TEST_CASE("decimation keeps a tone below the new Nyquist", "[dsp]") {
  const Audio a = tone(48000, 2.0, 200.0);
  const Audio d = decimate(a, 8000);
  CHECK(d.sr >= 8000);
  CHECK(d.duration() == Approx(a.duration()).margin(0.01));
  double rms = 0;
  for (size_t i = d.x.size() / 4; i < 3 * d.x.size() / 4; ++i)
    rms += d.x[i] * d.x[i];
  rms = std::sqrt(rms / (d.x.size() / 2));
  CHECK(rms == Approx(0.5 / std::sqrt(2.0)).epsilon(0.05));
}

TEST_CASE("16-bit WAV round-trips", "[dsp]") {
  const auto path =
      std::filesystem::temp_directory_path() / "audio2loc_test_roundtrip.wav";
  const Audio a = tone(8000, 1.0, 440.0, 0.25);
  write_wav16(path.string(), a);
  const Audio b = load_audio(path.string(), 8000);
  std::filesystem::remove(path);
  REQUIRE(b.sr == Approx(8000));
  REQUIRE(b.x.size() == a.x.size());
  for (size_t i = 0; i < a.x.size(); i += 97)
    CHECK(b.x[i] == Approx(a.x[i]).margin(2e-4));
}

TEST_CASE("the log-frequency image puts harmonics a constant row apart",
          "[dsp]") {
  Spectrogram s;
  s.nframes = 1;
  s.nbins = 200;
  s.df = 2.0;
  s.v.assign(200, 0.f);
  LogImage li = to_log_image(s, 50, 390, 0.01);
  const double r1 = li.row_of(100), r2 = li.row_of(200), r3 = li.row_of(300);
  CHECK(r2 - r1 == Approx(std::log(2.0) / 0.01));
  CHECK(r3 - r1 == Approx(std::log(3.0) / 0.01));
  CHECK(li.freq(r1) == Approx(100.0));
}
