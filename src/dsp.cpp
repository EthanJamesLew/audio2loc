#include "dsp.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "geo.h"
#include "gpu.h"
#include "util.h"

namespace pl {

// ---------------------------------------------------------------- WAV ---
static uint32_t rd32(const unsigned char* p) { return p[0] | (p[1] << 8) | (p[2] << 16) | ((uint32_t)p[3] << 24); }
static uint16_t rd16(const unsigned char* p) { return (uint16_t)(p[0] | (p[1] << 8)); }

static bool read_wav(const std::string& path, Audio& out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  std::vector<unsigned char> buf((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  if (buf.size() < 12 || std::memcmp(buf.data(), "RIFF", 4) != 0 || std::memcmp(buf.data() + 8, "WAVE", 4) != 0)
    return false;
  size_t p = 12;
  int fmt_tag = 0, ch = 0, bits = 0;
  uint32_t sr = 0;
  const unsigned char* data = nullptr;
  size_t data_len = 0;
  while (p + 8 <= buf.size()) {
    uint32_t len = rd32(&buf[p + 4]);
    const unsigned char* c = &buf[p + 8];
    size_t avail = std::min<size_t>(len, buf.size() - p - 8);
    if (!std::memcmp(&buf[p], "fmt ", 4) && avail >= 16) {
      fmt_tag = rd16(c);
      ch = rd16(c + 2);
      sr = rd32(c + 4);
      bits = rd16(c + 14);
      if (fmt_tag == 0xFFFE && avail >= 26) fmt_tag = rd16(c + 24);  // WAVE_FORMAT_EXTENSIBLE subformat
    } else if (!std::memcmp(&buf[p], "data", 4)) {
      data = c;
      data_len = avail;
    }
    p += 8 + len + (len & 1);
  }
  if (!data || ch <= 0 || sr == 0) die("%s: malformed WAV", path.c_str());
  int bps = bits / 8;
  size_t n = data_len / (bps * ch);
  out.sr = sr;
  out.x.assign(n, 0.f);
  for (size_t i = 0; i < n; ++i) {
    double acc = 0;
    for (int k = 0; k < ch; ++k) {
      const unsigned char* s = data + (i * ch + k) * bps;
      double v = 0;
      if (fmt_tag == 3 && bits == 32) {
        float fv;
        std::memcpy(&fv, s, 4);
        v = fv;
      } else if (fmt_tag == 3 && bits == 64) {
        double dv;
        std::memcpy(&dv, s, 8);
        v = dv;
      } else if (fmt_tag == 1 && bits == 8) {
        v = (s[0] - 128) / 128.0;
      } else if (fmt_tag == 1 && bits == 16) {
        v = (int16_t)rd16(s) / 32768.0;
      } else if (fmt_tag == 1 && bits == 24) {
        int32_t x = (s[0] << 8) | (s[1] << 16) | ((int32_t)s[2] << 24);
        v = (x >> 8) / 8388608.0;
      } else if (fmt_tag == 1 && bits == 32) {
        v = (int32_t)rd32(s) / 2147483648.0;
      } else {
        die("%s: unsupported WAV encoding (format %d, %d bits)", path.c_str(), fmt_tag, bits);
      }
      acc += v;
    }
    out.x[i] = (float)(acc / ch);
  }
  return true;
}

void write_wav16(const std::string& path, const Audio& a) {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) die("cannot write %s", path.c_str());
  uint32_t n = (uint32_t)a.x.size(), sr = (uint32_t)std::lround(a.sr);
  auto w32 = [&](uint32_t v) { std::fwrite(&v, 4, 1, f); };
  auto w16 = [&](uint16_t v) { std::fwrite(&v, 2, 1, f); };
  std::fwrite("RIFF", 1, 4, f);
  w32(36 + n * 2);
  std::fwrite("WAVEfmt ", 1, 8, f);
  w32(16);
  w16(1);
  w16(1);
  w32(sr);
  w32(sr * 2);
  w16(2);
  w16(16);
  std::fwrite("data", 1, 4, f);
  w32(n * 2);
  std::vector<int16_t> pcm(n);
  for (uint32_t i = 0; i < n; ++i) pcm[i] = (int16_t)std::lround(std::clamp(a.x[i], -1.f, 1.f) * 32767.f);
  std::fwrite(pcm.data(), 2, n, f);
  std::fclose(f);
}

Audio decimate(const Audio& a, double target_sr) {
  int D = (int)std::floor(a.sr / target_sr);
  if (D <= 1) return a;
  // Blackman-windowed sinc low-pass at 0.45 * output Nyquist... keep 90% of the new band.
  int half = 8 * D;
  double fc = 0.5 / D * 0.9;
  std::vector<float> h(2 * half + 1);
  double sum = 0;
  for (int k = -half; k <= half; ++k) {
    double x = k;
    double s = k == 0 ? 2 * fc : std::sin(2 * kPi * fc * x) / (kPi * x);
    double w = 0.42 + 0.5 * std::cos(kPi * x / half) + 0.08 * std::cos(2 * kPi * x / half);
    h[k + half] = (float)(s * w);
    sum += s * w;
  }
  for (auto& v : h) v = (float)(v / sum);
  Audio o;
  o.sr = a.sr / D;
  size_t n = a.x.size() / D;
  o.x.resize(n);
  const float* x = a.x.data();
  long N = (long)a.x.size();
  parallel_for(
      n,
      [&](size_t i, int) {
        long c = (long)i * D;
        double acc = 0;
        for (int k = -half; k <= half; ++k) {
          long j = c - k;
          if (j >= 0 && j < N) acc += h[k + half] * x[j];
        }
        o.x[i] = (float)acc;
      },
      4096);
  return o;
}

Audio load_audio(const std::string& path, double target_sr, const std::string& ffmpeg) {
  Audio a;
  if (!file_exists(path)) die("audio file not found: %s", path.c_str());
  if (read_wav(path, a)) {
    log("audio: %s  %.0f Hz  %.1f s (native WAV)", path.c_str(), a.sr, a.duration());
    return decimate(a, target_sr);
  }
  // Decode anything else with ffmpeg: mono float32 at a rate that is an integer multiple of target.
  int out_sr = (int)std::lround(target_sr);
  std::string esc;
  for (char c : path) {
    if (c == '\'') esc += "'\\''";
    else esc += c;
  }
  std::string cmd = ffmpeg + " -nostdin -v error -i '" + esc + "' -vn -ac 1 -ar " + std::to_string(out_sr) +
                    " -f f32le -acodec pcm_f32le - ";
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) die("failed to run ffmpeg (needed for non-WAV input). Install ffmpeg or pass --ffmpeg PATH");
  std::vector<float> buf(1 << 16);
  size_t got;
  while ((got = std::fread(buf.data(), sizeof(float), buf.size(), p)) > 0) a.x.insert(a.x.end(), buf.begin(), buf.begin() + got);
  int rc = pclose(p);
  if (rc != 0 || a.x.empty())
    die("ffmpeg failed to decode %s (exit %d). Is ffmpeg installed? Use --ffmpeg PATH or convert to WAV.",
        path.c_str(), rc);
  a.sr = out_sr;
  log("audio: %s  %d Hz  %.1f s (via ffmpeg)", path.c_str(), out_sr, a.duration());
  return a;
}

// ---------------------------------------------------------------- FFT ---
FFT::FFT(int n) : n_(n), rev_(n), tw_(n / 2) {
  if (n & (n - 1)) die("FFT size must be a power of two");
  int lg = 0;
  while ((1 << lg) < n) ++lg;
  for (int i = 0; i < n; ++i) {
    int r = 0;
    for (int b = 0; b < lg; ++b)
      if (i & (1 << b)) r |= 1 << (lg - 1 - b);
    rev_[i] = r;
  }
  for (int i = 0; i < n / 2; ++i) tw_[i] = std::polar(1.0f, (float)(-2.0 * kPi * i / n));
}

void FFT::forward(std::complex<float>* d) const {
  for (int i = 0; i < n_; ++i)
    if (i < rev_[i]) std::swap(d[i], d[rev_[i]]);
  for (int len = 2; len <= n_; len <<= 1) {
    int half = len >> 1, step = n_ / len;
    for (int i = 0; i < n_; i += len) {
      for (int j = 0; j < half; ++j) {
        std::complex<float> u = d[i + j], v = d[i + j + half] * tw_[j * step];
        d[i + j] = u + v;
        d[i + j + half] = u - v;
      }
    }
  }
}

// -------------------------------------------------------- spectrogram ---
Spectrogram compute_spectrogram(const Audio& a, double win_s, double hop_s, double t_begin, double t_end,
                                double fmin, double fmax, int zero_pad) {
  Spectrogram s;
  int wl = std::max(16, (int)std::lround(win_s * a.sr));
  int nfft = 1;
  while (nfft < wl * zero_pad) nfft <<= 1;
  s.win = wl / a.sr;
  s.hop = hop_s;
  s.df = a.sr / nfft;
  s.bin0 = std::max(0, (int)std::floor(fmin / s.df));
  int bin1 = std::min(nfft / 2, (int)std::ceil(fmax / s.df));
  s.nbins = std::max(1, bin1 - s.bin0 + 1);
  t_begin = std::max(t_begin, 0.0);
  t_end = std::min(t_end, a.duration());
  s.t0 = t_begin;
  s.nframes = std::max(0, (int)std::floor((t_end - t_begin) / hop_s) + 1);
  s.v.assign((size_t)s.nframes * s.nbins, 0.f);
  // GPU (cuFFT) for large transforms; identical windowing/normalisation
  if (gpu::enabled() && (double)s.nframes * nfft > 4e6 &&
      gpu::stft_db(a.x.data(), a.x.size(), a.sr, wl, nfft, s.t0, hop_s, s.nframes, s.bin0, s.nbins, s.v.data()))
    return s;
  std::vector<float> win(wl);
  double wsum = 0;
  for (int i = 0; i < wl; ++i) {
    win[i] = (float)(0.5 - 0.5 * std::cos(2 * kPi * (i + 0.5) / wl));
    wsum += win[i];
  }
  FFT fft(nfft);
  long N = (long)a.x.size();
  int nt = num_threads();
  std::vector<std::vector<std::complex<float>>> bufs(nt, std::vector<std::complex<float>>(nfft));
  double norm = 2.0 / wsum;
  parallel_for(
      (size_t)s.nframes,
      [&](size_t f, int tid) {
        auto& b = bufs[tid];
        long c = std::lround(s.time((int)f) * a.sr);
        long st = c - wl / 2;
        for (int i = 0; i < nfft; ++i) b[i] = 0;
        for (int i = 0; i < wl; ++i) {
          long j = st + i;
          if (j >= 0 && j < N) b[i] = a.x[j] * win[i];
        }
        fft.forward(b.data());
        float* row = &s.v[f * s.nbins];
        for (int k = 0; k < s.nbins; ++k) {
          double m = std::abs(b[s.bin0 + k]) * norm;
          row[k] = (float)(20.0 * std::log10(m + 1e-12));
        }
      },
      8);
  return s;
}

void whiten(Spectrogram& s, float clip) {
  parallel_for((size_t)s.nbins, [&](size_t b, int) {
    std::vector<float> col(s.nframes);
    for (int f = 0; f < s.nframes; ++f) col[f] = s.at(f, (int)b);
    if (col.empty()) return;
    auto mid = col.begin() + col.size() / 2;
    std::nth_element(col.begin(), mid, col.end());
    float med = *mid;
    for (auto& v : col) v = std::fabs(v - med);
    std::nth_element(col.begin(), mid, col.end());
    float mad = std::max(1e-3f, *mid * 1.4826f);
    for (int f = 0; f < s.nframes; ++f) s.at(f, (int)b) = std::clamp((s.at(f, (int)b) - med) / mad, 0.f, clip);
  });
}

LogImage to_log_image(const Spectrogram& w, double fmin, double fmax, double dln) {
  LogImage L;
  L.t0 = w.t0;
  L.hop = w.hop;
  L.nf = w.nframes;
  L.dln = dln;
  L.lnf0 = std::log(fmin);
  L.nr = std::max(1, (int)std::ceil((std::log(fmax) - L.lnf0) / dln) + 1);
  L.v.assign((size_t)L.nf * L.nr, 0.f);
  std::vector<int> b0(L.nr);
  std::vector<float> fr(L.nr);
  for (int r = 0; r < L.nr; ++r) {
    double x = L.freq(r) / w.df - w.bin0;
    int i = (int)std::floor(x);
    b0[r] = i;
    fr[r] = (float)(x - i);
  }
  parallel_for((size_t)L.nf, [&](size_t f, int) {
    float* o = &L.v[f * L.nr];
    for (int r = 0; r < L.nr; ++r) {
      int i = b0[r];
      if (i < 0 || i + 1 >= w.nbins) continue;
      o[r] = w.at((int)f, i) * (1 - fr[r]) + w.at((int)f, i + 1) * fr[r];
    }
  });
  return L;
}

LogImage LogImage::blurred(int rf, int rr) const {
  LogImage o = *this;
  std::vector<float> tmp(v.size());
  // along rows (frequency)
  for (int f = 0; f < nf; ++f) {
    const float* s = &v[(size_t)f * nr];
    float* d = &tmp[(size_t)f * nr];
    double acc = 0;
    int cnt = 0;
    for (int r = -rr; r < nr + rr; ++r) {
      int add = r + rr, sub = r - rr - 1;
      if (add >= 0 && add < nr) acc += s[add], ++cnt;
      if (sub >= 0 && sub < nr) acc -= s[sub], --cnt;
      if (r >= 0 && r < nr) d[r] = cnt ? (float)(acc / (2 * rr + 1)) : 0.f;
    }
  }
  // along frames (time); use max-of-window-like smoothing via mean
  for (int r = 0; r < nr; ++r) {
    double acc = 0;
    for (int f = -rf; f < nf + rf; ++f) {
      int add = f + rf, sub = f - rf - 1;
      if (add >= 0 && add < nf) acc += tmp[(size_t)add * nr + r];
      if (sub >= 0 && sub < nf) acc -= tmp[(size_t)sub * nr + r];
      if (f >= 0 && f < nf) o.v[(size_t)f * nr + r] = (float)(acc / (2 * rf + 1));
    }
  }
  return o;
}

}  // namespace pl
