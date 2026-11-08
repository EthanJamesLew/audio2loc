// Ground-truth scenario generator for validating the pipeline end to end.
//
// Produces: <prefix>.csv  synthetic ADS-B for a region/day (props + jets + decoy pairs),
//           <prefix>.wav  audio heard at a hidden observer: two prop flyovers `gap` seconds apart
//                         (one manoeuvring through its CPA), physically synthesised by evaluating each
//                         source's phase at the retarded emission time, 1/r spreading, in noise,
//           <prefix>_truth.txt  the hidden observer + the event boxes to pass to `locate`.
#include "simulate.h"

#include <fstream>
#include <random>

#include "adsb.h"
#include "doppler.h"
#include "dsp.h"
#include "geo.h"
#include "util.h"

namespace pl {
namespace {

struct Man {
  double t0, t1, turn_dps, accel, vz;
};
struct Truth {
  std::vector<double> t, lat, lon, alt, gs, trk, vz;
};

double hdg_rate(const std::vector<Man>& m, double t) {
  double r = 0;
  for (auto& x : m)
    if (t >= x.t0 && t < x.t1) r += x.turn_dps;
  return r;
}
double acc(const std::vector<Man>& m, double t) {
  double r = 0;
  for (auto& x : m)
    if (t >= x.t0 && t < x.t1) r += x.accel;
  return r;
}
double vzf(const std::vector<Man>& m, double t) {
  double r = 0;
  for (auto& x : m)
    if (t >= x.t0 && t < x.t1) r += x.vz;
  return r;
}

// Integrates kinematics from an anchor state forward to t_end and backward to t_start (1 s output).
Truth fly(double lat, double lon, double alt, double t_anchor, double hdg, double spd, const std::vector<Man>& man,
          double t_start, double t_end) {
  auto run = [&](double dir, double t_lim, Truth& out) {
    double la = lat, lo = lon, al = alt, h = hdg, v = spd;
    const double sub = 0.1;
    for (long step = 0;; ++step) {
      double t = t_anchor + dir * step * sub;
      if (step % 10 == 0) {
        out.t.push_back(t), out.lat.push_back(la), out.lon.push_back(lo), out.alt.push_back(al);
        out.gs.push_back(v), out.trk.push_back(std::fmod(std::fmod(h, 360) + 360, 360)), out.vz.push_back(vzf(man, t));
      }
      if ((dir > 0 && t >= t_lim) || (dir < 0 && t <= t_lim)) break;
      double tm = t + dir * sub / 2;
      double hm = h + dir * hdg_rate(man, tm) * sub / 2, vm = std::max(20.0, v + dir * acc(man, tm) * sub / 2);
      la += dir * vm * std::cos(hm * kDeg) * sub / m_per_deg_lat(la);
      lo += dir * vm * std::sin(hm * kDeg) * sub / m_per_deg_lon(la);
      al += dir * vzf(man, tm) * sub;
      h += dir * hdg_rate(man, tm) * sub;
      v = std::max(20.0, v + dir * acc(man, tm) * sub);
    }
  };
  Truth fwd, bwd, out;
  run(+1, t_end, fwd);
  run(-1, t_start, bwd);
  for (size_t i = bwd.t.size(); i-- > 1;) {
    out.t.push_back(bwd.t[i]), out.lat.push_back(bwd.lat[i]), out.lon.push_back(bwd.lon[i]);
    out.alt.push_back(bwd.alt[i]), out.gs.push_back(bwd.gs[i]), out.trk.push_back(bwd.trk[i]), out.vz.push_back(bwd.vz[i]);
  }
  for (size_t i = 0; i < fwd.t.size(); ++i) {
    out.t.push_back(fwd.t[i]), out.lat.push_back(fwd.lat[i]), out.lon.push_back(fwd.lon[i]);
    out.alt.push_back(fwd.alt[i]), out.gs.push_back(fwd.gs[i]), out.trk.push_back(fwd.trk[i]), out.vz.push_back(fwd.vz[i]);
  }
  return out;
}

Track to_track(const Truth& T) {
  Track tr;
  for (size_t i = 0; i < T.t.size(); ++i) {
    TrackPoint p;
    p.t = T.t[i];
    p.p = geodetic_to_ecef(T.lat[i], T.lon[i], T.alt[i]);
    p.v = enu_to_ecef_dir(T.lat[i], T.lon[i], T.gs[i] * std::sin(T.trk[i] * kDeg), T.gs[i] * std::cos(T.trk[i] * kDeg),
                          T.vz[i]);
    p.lat = (float)T.lat[i];
    p.lon = (float)T.lon[i];
    tr.pts.push_back(p);
  }
  tr.t0 = tr.pts.front().t;
  tr.t1 = tr.pts.back().t;
  return tr;
}

struct SimFlight {
  std::string icao, cs, type;
  Truth truth;
  Track track;
  double f0 = 100;
  bool jet = false;
};


// Physically synthesised audio at observer O: each source's harmonic phase is evaluated at the
// retarded emission time (which produces the exact Doppler shift), with 1/r spreading, plus noise.
Audio synth_audio(const std::vector<std::pair<const Track*, double>>& sources, const Vec3& O, double rec0,
                  double rec1, const SimOptions& o, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> N(0, 1);
  Audio a;
  a.sr = o.sr;
  size_t ns = (size_t)((rec1 - rec0) * o.sr);
  a.x.assign(ns, 0.f);
  std::vector<size_t> src;
  for (size_t i = 0; i < sources.size(); ++i) {
    const Track& tr = *sources[i].first;
    for (double t = std::max(tr.t0, rec0 - 60); t < std::min(tr.t1, rec1 + 60); t += 5) {
      Vec3 p, v;
      tr.eval(t, p, v);
      if ((p - O).norm() < 7000) {
        src.push_back(i);
        break;
      }
    }
  }
  log("simulate: %zu audible sources during the recording", src.size());
  const double c = o.c;
  const double harm[4] = {1.0, 0.55, 0.35, 0.2};
  for (size_t si : src) {
    const Track& track = *sources[si].first;
    double f0 = sources[si].second;
    double ph[4];
    for (auto& x : ph) x = U(rng) * 2 * kPi;
    size_t chunk = 1 << 16;
    parallel_for((ns + chunk - 1) / chunk, [&](size_t ci, int) {
      size_t a0 = ci * chunk, a1 = std::min(ns, a0 + chunk);
      int hint = -1;
      Vec3 p, v;
      double te = kNaN;
      for (size_t n = a0; n < a1; ++n) {
        double tr = rec0 + n / o.sr;
        if (!std::isfinite(te)) te = tr - 3;
        for (int k = 0; k < (n == a0 ? 10 : 2); ++k) {
          track.eval(te, p, v, &hint);
          te = tr - (p - O).norm() / c;
        }
        if (te < track.t0 || te > track.t1) continue;
        double r = (p - O).norm();
        double amp = o.amp_ref * 100.0 / std::max(r, 50.0);
        double tl = te - rec0;
        double phase0 = 2 * kPi * f0 * (tl + 0.002 * 40 / (2 * kPi) * std::sin(2 * kPi * tl / 40));  // rpm wobble
        double s = 0;
        for (int h = 0; h < 4; ++h) s += harm[h] * std::sin((h + 1) * phase0 + ph[h]);
        a.x[n] += (float)(amp * s);
      }
    });
  }
  for (auto& x : a.x) x += (float)(N(rng) * o.noise);
  float mx = 0;
  for (auto x : a.x) mx = std::max(mx, std::fabs(x));
  if (mx > 0.99f)
    for (auto& x : a.x) x *= 0.99f / mx;
  return a;
}

}  // namespace

int run_simulate(const SimOptions& o) {
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> N(0, 1);
  auto uni = [&](double a, double b) { return a + (b - a) * U(rng); };
  double T0 = parse_time(o.start), T1 = T0 + o.hours * 3600;
  double olat = o.lat + uni(-0.05, 0.05), olon = o.lon + uni(-0.05, 0.05);
  Vec3 O = geodetic_to_ecef(olat, olon, 0 + 2.0);

  std::vector<SimFlight> fl;
  int icao_seq = 0xa00000;
  auto mk_id = [&] {
    char b[16];
    std::snprintf(b, sizeof b, "%06x", icao_seq++);
    return std::string(b);
  };
  const char* props[] = {"C172", "P28A", "SR22", "BE36", "PC12", "DH8D", "C208", "BE20"};
  const char* jets[] = {"A320", "B738", "E75L", "CRJ9", "C56X", "B712"};

  // --- the two real events
  double tA = T0 + o.hours * 3600 * 0.45, tB = tA + o.gap;
  auto through = [&](double t_anchor, double hdg, double spd, double alt, double off, std::vector<Man> man,
                     const char* type, double f0) {
    // anchor = point offset perpendicular to heading from O
    double ph = (hdg + 90) * kDeg;
    double la = olat + off * std::cos(ph) / m_per_deg_lat(olat), lo = olon + off * std::sin(ph) / m_per_deg_lon(olat);
    SimFlight f;
    f.icao = mk_id();
    f.cs = "N" + std::to_string(100 + fl.size()) + "PL";
    f.type = type;
    f.f0 = f0;
    f.truth = fly(la, lo, alt, t_anchor, hdg, spd, man, t_anchor - 1500, t_anchor + 1500);
    f.track = to_track(f.truth);
    fl.push_back(std::move(f));
  };
  // A: ~100 m/s turboprop, banks through 40 degrees just after CPA and slows - not an ideal Doppler curve.
  through(tA, uni(0, 360), 100, 750, 520, {{tA + 4, tA + 24, 2.0, -0.4, 0}, {tA - 400, tA - 380, -1.5, 0, 0}}, "PC12", 112);
  // B: ~70 m/s piston single, gentle climbing turn before CPA.
  through(tB, uni(0, 360), 70, 480, 950, {{tB - 45, tB - 10, -1.2, 0.2, 2.0}}, "C182", 88);

  // --- decoy pairs: same speeds/distances/time offset elsewhere, straight and level
  for (int i = 0; i < o.n_decoy; ++i) {
    double qlat, qlon;
    do {
      qlat = o.lat + uni(-o.region / 2, o.region / 2);
      qlon = o.lon + uni(-o.region / 2, o.region / 2);
    } while (haversine_m(qlat, qlon, olat, olon) < 15000);
    double tq = uni(T0 + 1800, T1 - 1800 - o.gap);
    for (int k = 0; k < 2; ++k) {
      double spd = k ? 70 + uni(-2, 2) : 100 + uni(-3, 3);
      double alt = k ? 480 + uni(-60, 60) : 750 + uni(-80, 80);
      double off = k ? 950 + uni(-120, 120) : 520 + uni(-120, 120);
      double hdg = uni(0, 360), ph = (hdg + 90) * kDeg;
      double t = tq + (k ? o.gap + uni(-2, 2) : 0);
      double la = qlat + off * std::cos(ph) / m_per_deg_lat(qlat), lo = qlon + off * std::sin(ph) / m_per_deg_lon(qlat);
      SimFlight f;
      f.icao = mk_id();
      f.cs = "D" + std::to_string(i) + (k ? "B" : "A");
      f.type = k ? "C182" : "PC12";
      f.truth = fly(la, lo, alt, t, hdg, spd, {}, t - 1200, t + 1200);
      f.track = to_track(f.truth);
      fl.push_back(std::move(f));
    }
  }

  // --- background traffic (kept away from the observer during the recording)
  double rec0 = tA - o.lead, rec1 = tB + o.lead;
  auto near_obs = [&](const Truth& T) {
    for (size_t i = 0; i < T.t.size(); i += 5)
      if (T.t[i] > rec0 - 300 && T.t[i] < rec1 + 300 && haversine_m(T.lat[i], T.lon[i], olat, olon) < 6000) return true;
    return false;
  };
  for (int i = 0; i < o.n_prop + o.n_jet; ++i) {
    bool jet = i >= o.n_prop;
    for (int tries = 0; tries < 50; ++tries) {
      double dur = uni(1200, 4800);
      double ts = uni(T0, T1 - dur);
      double la = o.lat + uni(-o.region / 2, o.region / 2), lo = o.lon + uni(-o.region / 2, o.region / 2);
      double spd = jet ? uni(120, 240) : uni(45, 115), alt = jet ? uni(1500, 7000) : uni(300, 3000);
      std::vector<Man> man;
      for (int m = 0; m < 3; ++m) {
        double tm = uni(ts, ts + dur);
        man.push_back({tm, tm + uni(10, 40), uni(-3, 3), uni(-0.3, 0.3), uni(-3, 3)});
      }
      Truth T = fly(la, lo, alt, ts, uni(0, 360), spd, man, ts, ts + dur);
      if (near_obs(T)) continue;
      SimFlight f;
      f.icao = mk_id();
      f.cs = (jet ? "JET" : "N") + std::to_string(i);
      f.type = jet ? jets[i % 6] : props[i % 8];
      f.jet = jet;
      f.truth = std::move(T);
      f.track = to_track(f.truth);
      fl.push_back(std::move(f));
      break;
    }
  }

  // --- ADS-B CSV (irregular sampling, quantisation, position noise, occasional dropouts)
  {
    std::ofstream c(o.prefix + ".csv");
    c << "time,icao24,callsign,typecode,lat,lon,alt_ft,gs_kt,track,vrate_fpm\n";
    size_t rows = 0;
    for (auto& f : fl) {
      double t = f.truth.t.front() + uni(0, 2);
      int hint = -1;
      while (t < f.truth.t.back()) {
        Vec3 p, v;
        f.track.eval(t, p, v, &hint);
        double la, lo, h;
        ecef_to_geodetic(p, la, lo, h);
        double e, n, u;
        ecef_dir_to_enu(la, lo, v, e, n, u);
        la += N(rng) * 4 / m_per_deg_lat(la);
        lo += N(rng) * 4 / m_per_deg_lon(la);
        double alt_ft = std::round(h / kFt / 25) * 25;
        double gs = std::round(std::hypot(e, n) / kKt);
        double trk = std::fmod(std::atan2(e, n) / kDeg + 360, 360);
        double vr = std::round(u / kFt * 60 / 64) * 64;
        char b[256];
        std::snprintf(b, sizeof b, "%.2f,%s,%s,%s,%.6f,%.6f,%.0f,%.0f,%.1f,%.0f\n", t, f.icao.c_str(), f.cs.c_str(),
                      f.type.c_str(), la, lo, alt_ft, gs, trk, vr);
        c << b;
        ++rows;
        t += U(rng) < 0.01 ? uni(8, 25) : uni(0.5, 3.0);
      }
    }
    log("simulate: %zu flights, %zu ADS-B rows -> %s.csv", fl.size(), rows, o.prefix.c_str());
  }

  // --- audio at the observer
  std::vector<std::pair<const Track*, double>> srcs;
  for (auto& f : fl) srcs.push_back({&f.track, f.f0});
  Audio a = synth_audio(srcs, O, rec0, rec1, o, rng);
  write_wav16(o.prefix + ".wav", a);
  const double c = o.c;

  // --- truth + event boxes
  std::ofstream tf(o.prefix + "_truth.txt");
  char b[1024];
  std::snprintf(b, sizeof b, "observer_lat %.6f\nobserver_lon %.6f\nrecording_start %s\n", olat, olon,
                fmt_time(rec0).c_str());
  tf << b;
  std::string evargs;
  for (int k = 0; k < 2; ++k) {
    double tcpa, d, sp;
    double tg = k ? tB : tA;
    // brute-force CPA on the truth track
    double best = 1e30;
    tcpa = tg;
    for (double t = tg - 200; t < tg + 200; t += 0.05) {
      Vec3 p, v;
      fl[k].track.eval(t, p, v);
      double dd = (p - O).norm();
      if (dd < best) best = dd, tcpa = t;
    }
    d = best;
    Vec3 p, v;
    fl[k].track.eval(tcpa, p, v);
    sp = v.norm();
    double tc_audio = tcpa + d / c - rec0;
    std::snprintf(b, sizeof b, "event%d flight %s %s  cpa_emit %s  slant %.0f m  speed %.1f m/s  audio_tc %.2f s  f0 %.1f Hz\n",
                  k + 1, fl[k].icao.c_str(), fl[k].type.c_str(), fmt_time(tcpa).c_str(), d, sp, tc_audio, fl[k].f0);
    tf << b;
    std::snprintf(b, sizeof b, " --event t0=%.0f,t1=%.0f,fmin=%.0f,fmax=%.0f,mult=1:2:3", tc_audio - o.box,
                  tc_audio + o.box, fl[k].f0 * 0.72, fl[k].f0 * 1.42);
    evargs += b;
  }
  std::snprintf(b, sizeof b, "suggested: planeloc locate --audio %s.wav --adsb %s.csv%s --start %s --end %s\n",
                o.prefix.c_str(), o.prefix.c_str(), evargs.c_str(), fmt_time(T0).c_str(), fmt_time(T1).c_str());
  tf << b;
  tf.close();
  std::ifstream rd(o.prefix + "_truth.txt");
  std::string line;
  while (std::getline(rd, line)) std::printf("%s\n", line.c_str());
  return 0;
}

}  // namespace pl

namespace pl {

// Scenario on REAL traffic: pick a hidden observer under a real (non-jet) flight, find the next real
// non-jet flight to pass close by, require silence in between, synthesise the audio from those real
// trajectories. Searching it with `locate` exercises the full real-world data volume.
int run_simulate_real(const SimOptions& o, const std::string& adsb) {
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0, 1);
  LoadOptions lo;
  std::vector<Flight> flights;
  load_adsb(adsb, lo, flights);
  auto tracks = build_tracks(flights, 90);
  std::vector<size_t> elig;
  for (size_t i = 0; i < tracks.size(); ++i)
    if (tracks[i].cls != ACClass::Jet && tracks[i].t1 - tracks[i].t0 > 1200) elig.push_back(i);
  log("simulate-real: %zu tracks, %zu eligible", tracks.size(), elig.size());
  struct Pass {
    size_t tr;
    double trx, d;
  };
  for (int attempt = 0; attempt < 400; ++attempt) {
    const Track& T = tracks[elig[(size_t)(U(rng) * elig.size())]];
    double t = T.t0 + 300 + U(rng) * (T.t1 - T.t0 - 600);
    Vec3 p, v;
    T.eval(t, p, v);
    double lat, lon, h;
    ecef_to_geodetic(p, lat, lon, h);
    double spd = v.norm();
    if (h < 200 || h > 1500 || spd < 40 || spd > 130) continue;
    double e, n, u;
    ecef_dir_to_enu(lat, lon, v, e, n, u);
    double hd = std::atan2(e, n) + (U(rng) < 0.5 ? 1 : -1) * kPi / 2, off = 250 + U(rng) * 600;
    double olat = lat + off * std::cos(hd) / m_per_deg_lat(lat), olon = lon + off * std::sin(hd) / m_per_deg_lon(lat);
    Vec3 O = geodetic_to_ecef(olat, olon, 2.0);
    // all passes within 5 km of O around this time
    std::vector<Pass> passes;
    for (size_t i = 0; i < tracks.size(); ++i) {
      const Track& R = tracks[i];
      if (R.t1 < t - 1800 || R.t0 > t + 2700 || R.maxlat < olat - 0.1 || R.minlat > olat + 0.1 ||
          R.maxlon < olon - 0.15 || R.minlon > olon + 0.15)
        continue;
      double prev = 1e30, prev2 = 1e30;
      int hint = -1;
      for (double tt = std::max(R.t0, t - 1800); tt <= std::min(R.t1, t + 2700); tt += 2) {
        R.eval(tt, p, v, &hint);
        double d = (p - O).norm();
        if (prev < prev2 && prev <= d && prev < 5000) {
          double tc, dc, sp;
          if (track_cpa(R, O, tt - 6, tt + 2, tc, dc, sp)) passes.push_back({i, tc + dc / o.c, dc});
        }
        prev2 = prev;
        prev = d;
      }
    }
    std::sort(passes.begin(), passes.end(), [](const Pass& a, const Pass& b) { return a.trx < b.trx; });
    // event 1 = the pass of T closest to t; event 2 = next pass within 2 km (non-jet), quiet in between
    int i1 = -1;
    for (size_t k = 0; k < passes.size(); ++k)
      if (&tracks[passes[k].tr] == &T && std::fabs(passes[k].trx - t) < 120 && passes[k].d < o.max_slant) i1 = (int)k;
    if (i1 < 0) continue;
    int i2 = -1;
    bool clean = true;
    for (size_t k = i1 + 1; k < passes.size(); ++k) {
      const Pass& q = passes[k];
      bool jet = tracks[q.tr].cls == ACClass::Jet;
      if (q.d < (jet ? 4500 : 2500)) {
        if (!jet && q.d < o.max_slant && q.trx - passes[i1].trx > 240) i2 = (int)k;
        else clean = false;
        break;
      }
    }
    for (int k = i1 - 1; k >= 0; --k)
      if (passes[i1].trx - passes[k].trx < 300 && passes[k].d < 3000) clean = false;
    if (i2 < 0 || !clean || passes[i2].trx - passes[i1].trx > 2400) continue;
    if (i2 + 1 < (int)passes.size() && passes[i2 + 1].trx - passes[i2].trx < 300 && passes[i2 + 1].d < 3000) continue;
    double T1 = passes[i1].trx, T2 = passes[i2].trx;
    double rec0 = T1 - 240 - U(rng) * 60, rec1 = T2 + 240;
    std::vector<std::pair<const Track*, double>> srcs;
    double f1 = 85 + U(rng) * 45, f2 = 75 + U(rng) * 45;
    for (auto& q : passes) {
      double f = q.tr == passes[i1].tr ? f1 : q.tr == passes[i2].tr ? f2 : 70 + U(rng) * 70;
      srcs.push_back({&tracks[q.tr], f});
    }
    Audio a = synth_audio(srcs, O, rec0, rec1, o, rng);
    write_wav16(o.prefix + ".wav", a);
    std::ofstream tf(o.prefix + "_truth.txt");
    char b[1024];
    std::snprintf(b, sizeof b, "observer_lat %.6f\nobserver_lon %.6f\nrecording_start %s\n", olat, olon,
                  fmt_time(rec0).c_str());
    tf << b;
    std::string ev;
    int k = 0;
    for (int idx : {i1, i2}) {
      const Pass& q = passes[idx];
      double fz = idx == i1 ? f1 : f2;
      std::snprintf(b, sizeof b, "event%d flight %s  slant %.0f m  heard %s  audio_tc %.2f s  f0 %.1f Hz\n", ++k,
                    track_label(flights, tracks[q.tr]).c_str(), q.d, fmt_time(q.trx).c_str(), q.trx - rec0, fz);
      tf << b;
      std::snprintf(b, sizeof b, " --event t0=%.0f,t1=%.0f,fmin=%.0f,fmax=%.0f,mult=1:2:3", q.trx - rec0 - o.box,
                    q.trx - rec0 + o.box, fz * 0.72, fz * 1.42);
      ev += b;
    }
    std::string day = fmt_time(std::floor(T1 / 86400) * 86400).substr(0, 10);
    std::snprintf(b, sizeof b, "suggested: planeloc locate --audio %s.wav --adsb %s%s --start %s --end %sT23:59:59Z\n",
                  o.prefix.c_str(), adsb.c_str(), ev.c_str(), day.c_str(), day.c_str());
    tf << b;
    tf.close();
    std::ifstream rd(o.prefix + "_truth.txt");
    std::string line;
    while (std::getline(rd, line)) std::printf("%s\n", line.c_str());
    return 0;
  }
  log("simulate-real: no suitable observer found");
  return 1;
}

}  // namespace pl
