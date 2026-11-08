// planeloc - geolocate a recording from the Doppler curves of aircraft flyovers + ADS-B history.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <tuple>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "adsb.h"
#include "dem.h"
#include "doppler.h"
#include "dsp.h"
#include "export.h"
#include "search.h"
#include "simulate.h"
#include "gpu.h"
#include "mcp.h"
#include "update.h"
#include "util.h"

using namespace pl;

static const char* kUsage = R"(planeloc - locate a recording from aircraft Doppler curves + ADS-B history

usage: planeloc <command> [options]

commands
  spectrogram  render a spectrogram PNG to find Doppler curves (zoom with --t0/--t1/--fmin/--fmax)
  measure      fit the boxed Doppler curve(s): speed, CPA slant distance, CPA time, f0 (+ overlay PNG)
  locate       full geolocation search -> ranked locations, GeoJSON layers, QGIS project
  ingest       convert ADS-B (CSV / readsb-adsb.lol trace dirs) into a fast binary cache (.plb)
  update       download/refresh the per-day ADS-B database from adsb.lol (global, streamed)
  mcp          run as a Model Context Protocol server on stdio (for AI agents)
  simulate     generate a ground-truth test scenario (ADS-B CSV + WAV + truth)

common
  --audio FILE          WAV natively; video/other audio via ffmpeg (--ffmpeg PATH)
  --event SPEC          one Doppler curve; repeat for each. SPEC is comma separated key=value:
                          t0,t1        audio time box (s)        fmin,fmax  frequency box (Hz) around ONE tone
                          mult=1:2:3   harmonics of that tone to score          side=both|before|after (fit)
                          score=brightness|ridge (ridge: trace-based, for very faint curves)
                          trace=FILE   hand-traced curve, CSV "time_s,freq_hz" (implies score=ridge)
                          overrides: v=MIN:MAX d=MIN:MAX tc=SEC tc_tol=SEC f0=MIN:MAX
  --win S --hop S       spectrogram window/hop (default 1.5 / 0.1)
  --temp C              air temperature for speed of sound (default 20 C)
  -j N                  threads          -v / -q   verbose / quiet
  --no-gpu              force the CPU path (CUDA builds use the GPU automatically)   --gpu-device N

locate
  --adsb PATH           ADS-B source (repeatable): .csv[.gz] | .plb | trace .json[.gz] | directory
  --db DIR              use the per-day database (adsb_YYYY-MM-DD.plb) for the days in the time window
  --auto-update         download missing days into --db first
  --start ISO --end ISO approximate date/time window to search (e.g. 2025-09-08 .. 2025-09-09)
  --rec-start ISO       approximate absolute time of audio t=0 (optional; tightens the search a lot)
  --rec-start-tol S     uncertainty of --rec-start (default 300 s)
  --bbox W,S,E,N        spatial bounds (deg) (default: everything in the data)
  --dem DIR             directory of SRTM-style .hgt tiles (N40W074.hgt); default flat --elev 0
  --elev M              ground elevation when no DEM tile           --obs-height M (default 2)
  --types FILE          extra type->class CSV (typecode,jet|turboprop|piston|heli)
  --include-jets        allow jets as event sources
  --grid M              coarse grid (default 200)  --fine M (default 100)
  --v-margin F --d-margin F --tc-margin S --offset-slack S   loosen measured ranges
  --no-silence          disable silence filter      --quiet A:B  extra audio interval with no aircraft (repeatable)
  --silence-radius M    (default 2000)  --silence-radius-jet M (4000)  --silence-guard S (60)
  --dt-sweep S          event-1 time alignment sweep (default 3)
  --refine-step M --refine-radius M --refine-top N   dense heatmap around winners (25 / 1500 / 1)
  --top N --overlays N  ranked rows / overlay PNGs (25 / 3)
  --max-alt M           ignore ADS-B above this altitude MSL (default 7000)
  --mem-limit GB        hard memory cap (default: 50% of available RAM; env PLANELOC_MEM_LIMIT). The search
                        works tile by tile and splits tiles to stay under it; a watchdog exits cleanly
                        (code 3) rather than letting the machine run out of memory
  --tile-deg D          initial spatial tile size for the candidate stages (default 2)
  --geoid N             geoid undulation (m) for geometric altitudes (CONUS ~ -30)
  --out DIR             output directory (default planeloc_out)
  --qgz FILE            also write a .qgz through the native QGIS API (needs -DPLANELOC_WITH_QGIS=ON)

ingest:   --adsb PATH... -o cache.plb [--start --end --bbox --max-alt]
update:   --db DIR (default $PLANELOC_DB or ./planeloc_db)  --date YYYY-MM-DD (repeatable) | --from D --to D |
          --latest | --list   [--bbox W,S,E,N --max-alt M --force]   (set GITHUB_TOKEN to lift API limits)
mcp:      --work DIR (default $PLANELOC_WORK or ./planeloc_work)  --db DIR
simulate: --prefix P --seed N --start ISO --hours H --lat --lon --region DEG --decoys N --props N --jets N
          --noise X --gap S
)";

struct Args {
  std::string cmd;
  std::multimap<std::string, std::string> kv;
  bool has(const std::string& k) const { return kv.count(k) > 0; }
  std::string get(const std::string& k, const std::string& def = "") const {
    auto it = kv.find(k);
    return it == kv.end() ? def : it->second;
  }
  double num(const std::string& k, double def) const {
    if (!has(k)) return def;
    double v;
    if (!parse_double(get(k), v)) die("--%s expects a number, got '%s'", k.c_str(), get(k).c_str());
    return v;
  }
  std::vector<std::string> all(const std::string& k) const {
    std::vector<std::string> o;
    auto r = kv.equal_range(k);
    for (auto it = r.first; it != r.second; ++it) o.push_back(it->second);
    return o;
  }
};

static Args parse_args(int argc, char** argv) {
  static const char* flags[] = {"include-jets", "no-silence", "no-osm", "v", "q", "help", "h", "auto-update", "latest", "force", "list", "no-gpu"};
  Args a;
  if (argc < 2) {
    std::fputs(kUsage, stdout);
    std::exit(0);
  }
  a.cmd = argv[1];
  for (int i = 2; i < argc; ++i) {
    std::string s = argv[i];
    if (s.size() < 2 || s[0] != '-') die("unexpected argument '%s'", s.c_str());
    std::string k = s.substr(s[1] == '-' ? 2 : 1), v;
    size_t eq = k.find('=');
    if (eq != std::string::npos && k != "event") {
      v = k.substr(eq + 1);
      k = k.substr(0, eq);
    } else {
      bool is_flag = false;
      for (auto f : flags) is_flag |= k == f;
      if (!is_flag) {
        if (i + 1 >= argc) die("--%s needs a value", k.c_str());
        v = argv[++i];
      }
    }
    a.kv.emplace(k, v);
  }
  return a;
}

static std::pair<double, double> range2(const std::string& s) {
  auto p = split(s, ':');
  double a, b;
  if (p.size() != 2 || !parse_double(p[0], a) || !parse_double(p[1], b)) die("bad range '%s' (want MIN:MAX)", s.c_str());
  return {a, b};
}

static EventSpec parse_event(const std::string& spec) {
  EventSpec e;
  bool have[4] = {};
  for (auto& kvs : split(spec, ',')) {
    auto eq = kvs.find('=');
    if (eq == std::string::npos) die("bad --event item '%s'", kvs.c_str());
    std::string k = lower(trim(kvs.substr(0, eq))), v = trim(kvs.substr(eq + 1));
    double x = 0;
    auto num = [&] {
      if (!parse_double(v, x)) die("--event %s: not a number '%s'", k.c_str(), v.c_str());
      return x;
    };
    if (k == "t0") e.t0 = num(), have[0] = true;
    else if (k == "t1") e.t1 = num(), have[1] = true;
    else if (k == "fmin") e.fmin = num(), have[2] = true;
    else if (k == "fmax") e.fmax = num(), have[3] = true;
    else if (k == "mult") {
      e.mult.clear();
      for (auto& m : split(v, ':')) {
        double y;
        if (!parse_double(m, y) || y <= 0) die("bad mult '%s'", v.c_str());
        e.mult.push_back(y);
      }
    } else if (k == "side") e.side = v == "before" ? -1 : v == "after" ? 1 : 0;
    else if (k == "score") e.ridge_score = v == "ridge";
    else if (k == "trace") e.trace = v, e.ridge_score = true;
    else if (k == "v") std::tie(e.vmin, e.vmax) = range2(v);
    else if (k == "d") std::tie(e.dmin, e.dmax) = range2(v);
    else if (k == "f0") std::tie(e.f0min, e.f0max) = range2(v);
    else if (k == "tc") e.tc = num();
    else if (k == "tc_tol") e.tc_tol = num();
    else die("unknown --event key '%s'", k.c_str());
  }
  if (!(have[0] && have[1] && have[2] && have[3])) die("--event needs t0, t1, fmin, fmax");
  if (e.t1 <= e.t0 || e.fmax <= e.fmin) die("--event box is empty");
  return e;
}

static SearchParams common_params(const Args& a) {
  SearchParams P;
  P.c = speed_of_sound(a.num("temp", 20));
  P.win = a.num("win", P.win);
  P.hop = a.num("hop", P.hop);
  P.target_sr = a.num("sr", P.target_sr);
  return P;
}

static int cmd_spectrogram(const Args& a) {
  SearchParams P = common_params(a);
  Audio au = load_audio(a.get("audio"), P.target_sr, a.get("ffmpeg", "ffmpeg"));
  double t0 = a.num("t0", 0), t1 = a.num("t1", au.duration()), f0 = a.num("fmin", 20), f1 = a.num("fmax", 600);
  Spectrogram s = compute_spectrogram(au, P.win, P.hop, t0, t1, f0, f1);
  whiten(s);
  int ys = (int)a.num("yscale", std::max(1.0, std::round(400.0 / std::max(1, s.nbins))));
  Image img = render_spectrogram(s, t0, t1, f0, f1, {}, ys);
  std::string out = a.get("o", a.get("out", "spectrogram.png"));
  write_png(out, img);
  std::printf("%s: %dx%d px  x: %.3f s/px from t=%.2f  y: %.3f Hz/px (top=%.1f Hz, bottom=%.1f Hz)\n", out.c_str(), img.w,
              img.h, s.hop, s.t0, s.df / ys, f1, f0);
  return 0;
}

static int cmd_measure(const Args& a) {
  SearchParams P = common_params(a);
  Audio au = load_audio(a.get("audio"), P.target_sr, a.get("ffmpeg", "ffmpeg"));
  std::vector<EventSpec> specs;
  for (auto& s : a.all("event")) specs.push_back(parse_event(s));
  if (specs.empty()) die("measure needs at least one --event");
  auto ev = prepare_events(au, specs, P);
  print_events(ev);
  std::string out = a.get("out", "planeloc_measure");
  std::filesystem::create_directories(out);
  for (size_t k = 0; k < ev.size(); ++k) {
    const Event& e = ev[k];
    Curve ridge{{}, {60, 255, 60}}, fit{{}, {80, 255, 255}};
    for (size_t i = 0; i < e.ridge.t.size(); ++i) ridge.tf.push_back({e.ridge.t[i], e.ridge.f[i]});
    for (double t = e.spec.t0; t <= e.spec.t1; t += P.hop)
      fit.tf.push_back({t, e.fit.f0 * std::exp(straight_lng(t, e.fit.v, e.fit.d, e.fit.tc, P.c))});
    Image img = render_spectrogram(e.white, e.spec.t0 - 5, e.spec.t1 + 5, e.spec.fmin * 0.85, e.spec.fmax * 1.15,
                                   {ridge, fit}, 2);
    std::string p = out + "/measure_event" + std::to_string(k + 1) + ".png";
    write_png(p, img);
    std::printf("wrote %s (green: traced ridge, cyan: straight-line fit)\n", p.c_str());
  }
  return 0;
}

static LoadOptions load_opts(const Args& a) {
  LoadOptions o;
  if (a.has("start")) o.tmin = parse_time(a.get("start"));
  if (a.has("end")) o.tmax = parse_time(a.get("end"));
  if (a.has("start") && !a.has("end")) o.tmax = o.tmin + 86400;
  if (a.has("bbox")) {
    auto p = split(a.get("bbox"), ',');
    if (p.size() != 4) die("--bbox W,S,E,N");
    double v[4];
    for (int i = 0; i < 4; ++i)
      if (!parse_double(p[i], v[i])) die("--bbox W,S,E,N");
    double m = 0.15;  // keep flights just outside the box: they can be heard inside it
    o.minlon = v[0] - m, o.minlat = v[1] - m, o.maxlon = v[2] + m, o.maxlat = v[3] + m;
  }
  o.max_alt_m = a.num("max-alt", o.max_alt_m);
  o.geoid_undulation = a.num("geoid", 0);
  o.csv_alt_unit = a.get("csv-alt-unit");
  o.csv_speed_unit = a.get("csv-speed-unit");
  return o;
}

static int cmd_ingest(const Args& a) {
  start_mem_watchdog(resolve_mem_limit(a.num("mem-limit", 0)));
  LoadOptions o = load_opts(a);
  o.min_dt = a.num("min-dt", 0);  // keep full resolution in caches unless asked
  std::vector<Flight> fl;
  for (auto& p : a.all("adsb")) load_adsb(p, o, fl);
  std::string out = a.get("o", a.get("out", "adsb.plb"));
  save_cache(out, fl);
  size_t ns = 0;
  for (auto& f : fl) ns += f.s.size();
  std::printf("wrote %s: %zu flight legs, %zu samples\n", out.c_str(), fl.size(), ns);
  return 0;
}

static int cmd_locate(const Args& a) {
  SearchParams P = common_params(a);
  if (a.has("types")) load_type_table(a.get("types"));
  std::vector<EventSpec> specs;
  for (auto& s : a.all("event")) specs.push_back(parse_event(s));
  if (specs.empty()) die("locate needs at least one --event");
  if (!a.has("audio")) die("locate needs --audio");
  if (!a.has("adsb") && !a.has("db")) die("locate needs --adsb or --db");
  P.v_margin = a.num("v-margin", P.v_margin);
  P.d_margin = a.num("d-margin", P.d_margin);
  P.tc_margin = a.num("tc-margin", P.tc_margin);
  P.f0_margin = a.num("f0-margin", P.f0_margin);
  P.offset_slack = a.num("offset-slack", P.offset_slack);
  P.coarse_step = a.num("grid", P.coarse_step);
  P.fine_step = a.num("fine", P.fine_step);
  P.sample_dt = a.num("sample-dt", P.sample_dt);
  P.obs_height = a.num("obs-height", P.obs_height);
  P.silence = !a.has("no-silence");
  P.silence_radius = a.num("silence-radius", P.silence_radius);
  P.silence_radius_jet = a.num("silence-radius-jet", P.silence_radius_jet);
  P.silence_guard = a.num("silence-guard", P.silence_guard);
  for (auto& q : a.all("quiet")) P.quiet.push_back(range2(q));
  P.include_jets = a.has("include-jets");
  P.dt_sweep = a.num("dt-sweep", P.dt_sweep);
  P.top = (int)a.num("top", P.top);
  P.overlays = (int)a.num("overlays", P.overlays);
  P.refine_top = (int)a.num("refine-top", P.refine_top);
  P.refine_step = a.num("refine-step", P.refine_step);
  P.refine_radius = a.num("refine-radius", P.refine_radius);
  P.max_candidates = (size_t)a.num("max-candidates", (double)P.max_candidates);
  P.tile_deg = a.num("tile-deg", P.tile_deg);
  P.mem_limit = resolve_mem_limit(a.num("mem-limit", 0));
  start_mem_watchdog(P.mem_limit);
  log("memory limit %.1f GB (override with --mem-limit GB or PLANELOC_MEM_LIMIT)", gb(P.mem_limit));
  P.outdir = a.get("out", P.outdir);
  P.qgis_native_out = a.get("qgz");
  LoadOptions lo = load_opts(a);
  P.tmin = lo.tmin;
  P.tmax = lo.tmax;
  if (a.has("rec-start")) {
    P.rec_start = parse_time(a.get("rec-start"));
    P.rec_start_tol = a.num("rec-start-tol", P.rec_start_tol);
    double span = P.rec_start_tol + 4 * 3600;  // recordings are rarely > a few hours
    lo.tmin = std::max(lo.tmin, P.rec_start - span);
    lo.tmax = std::min(lo.tmax, P.rec_start + span);
  }

  Audio au = load_audio(a.get("audio"), P.target_sr, a.get("ffmpeg", "ffmpeg"));
  auto events = prepare_events(au, specs, P);
  print_events(events);

  std::vector<std::string> sources = a.all("adsb");
  if (a.has("db")) {
    double t0 = std::isfinite(P.rec_start) ? P.rec_start - P.rec_start_tol - 600 : lo.tmin;
    double t1 = std::isfinite(P.rec_start) ? P.rec_start + P.rec_start_tol + au.duration() + 600 : lo.tmax;
    if (t0 < -1e299 || t1 > 1e299) die("--db needs a time window (--start/--end or --rec-start)");
    UpdateOptions uo;
    uo.db = a.get("db");
    uo.dates = dates_between(t0, t1);
    if (uo.dates.size() > 14) die("time window spans %zu days; narrow --start/--end", uo.dates.size());
    if (a.has("auto-update")) {
      uo.filter.max_alt_m = lo.max_alt_m;
      run_update(uo);
    }
    for (auto& d : uo.dates) {
      std::string p = db_day_path(uo.db, d);
      if (file_exists(p)) {
        sources.push_back(p);
        double w, s, e, n;
        if (db_day_bbox(p, w, s, e, n) && (w > -179.9 || e < 179.9 || s > -89.9 || n < 89.9) &&
            (lo.minlon < w || lo.maxlon > e || lo.minlat < s || lo.maxlat > n))
          log("db: WARNING %s only covers bbox %.2f,%.2f,%.2f,%.2f - the search area extends outside it", p.c_str(),
              w, s, e, n);
      }
      else log("db: missing %s (run `planeloc update --db %s --date %s` or pass --auto-update)", p.c_str(),
               uo.db.c_str(), d.c_str());
    }
  }
  if (sources.empty()) die("no ADS-B source: pass --adsb or --db");
  std::vector<Flight> flights;
  for (auto& p : sources) load_adsb(p, lo, flights);
  if (flights.empty()) die("no ADS-B data in the requested window");
  auto tracks = build_tracks(flights, a.num("max-gap", 90));
  // raw samples are no longer needed (labels only from here on)
  size_t raw = 0;
  for (auto& f : flights) raw += f.s.capacity() * sizeof(RawSample), std::vector<RawSample>().swap(f.s);
  log("tracks: %zu legs built; released %.1f GB of raw samples (rss %.1f GB)", tracks.size(), gb(raw),
      gb(current_rss_bytes()));

  Dem dem(a.get("dem"), a.num("elev", 0));
  if (a.has("dem")) {
    double mila = 90, mala = -90, milo = 180, malo = -180;
    for (auto& t : tracks)
      mila = std::min(mila, t.minlat), mala = std::max(mala, t.maxlat), milo = std::min(milo, t.minlon),
      malo = std::max(malo, t.maxlon);
    dem.preload(mila - 0.1, mala + 0.1, milo - 0.1, malo + 0.1);
    if (!dem.has_data()) log("dem: no tiles found in %s - using constant %.0f m", a.get("dem").c_str(), a.num("elev", 0));
  }
  return run_locate(au, flights, tracks, dem, events, P);
}

static std::string env_or(const char* k, const char* def) {
  const char* v = std::getenv(k);
  return v && *v ? v : def;
}

static int cmd_update(const Args& a) {
  UpdateOptions o;
  o.db = a.get("db", env_or("PLANELOC_DB", "planeloc_db"));
  if (a.has("list")) {
    print_db(o.db);
    return 0;
  }
  LoadOptions lo = load_opts(a);
  o.filter.minlat = lo.minlat, o.filter.maxlat = lo.maxlat, o.filter.minlon = lo.minlon, o.filter.maxlon = lo.maxlon;
  o.filter.max_alt_m = lo.max_alt_m;
  o.filter.geoid_undulation = lo.geoid_undulation;
  o.force = a.has("force");
  o.curl = a.get("curl", "curl");
  for (auto& d : a.all("date")) o.dates.push_back(utc_date(parse_time(d)));
  if (a.has("from")) {
    double t0 = parse_time(a.get("from")), t1 = a.has("to") ? parse_time(a.get("to")) : t0;
    for (auto& d : dates_between(t0, t1)) o.dates.push_back(d);
  }
  if (a.has("latest") || o.dates.empty()) {
    std::string d = latest_available_date(o.curl);
    if (d.empty()) die("could not find a recent adsb.lol release (network? GitHub rate limit? set GITHUB_TOKEN)");
    log("update: newest published day is %s", d.c_str());
    o.dates.push_back(d);
  }
  int fails = run_update(o);
  print_db(o.db);
  return fails ? 1 : 0;
}

static int cmd_mcp(const Args& a) {
  g_verbose = 0;  // stdout belongs to the protocol
  std::string self = std::filesystem::read_symlink("/proc/self/exe").string();
  return run_mcp_server(self, a.get("work", env_or("PLANELOC_WORK", "planeloc_work")),
                        a.get("db", env_or("PLANELOC_DB", "planeloc_db")));
}

static int cmd_simulate(const Args& a) {
  SimOptions o;
  o.prefix = a.get("prefix", o.prefix);
  o.seed = (uint64_t)a.num("seed", (double)o.seed);
  o.start = a.get("start", o.start);
  o.hours = a.num("hours", o.hours);
  o.lat = a.num("lat", o.lat);
  o.lon = a.num("lon", o.lon);
  o.region = a.num("region", o.region);
  o.n_decoy = (int)a.num("decoys", o.n_decoy);
  o.n_prop = (int)a.num("props", o.n_prop);
  o.n_jet = (int)a.num("jets", o.n_jet);
  o.noise = a.num("noise", o.noise);
  o.amp_ref = a.num("amp", o.amp_ref);
  o.gap = a.num("gap", o.gap);
  o.max_slant = a.num("max-slant", o.max_slant);
  o.c = speed_of_sound(a.num("temp", 20));
  if (a.has("adsb")) {
    start_mem_watchdog(resolve_mem_limit(a.num("mem-limit", 0)));
    return run_simulate_real(o, a.get("adsb"));
  }
  return run_simulate(o);
}

int main(int argc, char** argv) {
  Args a = parse_args(argc, argv);
  if (a.has("v")) g_verbose = 2;
  if (a.has("q")) g_verbose = 0;
  if (a.has("j")) g_threads = (int)a.num("j", 0);
  if (a.cmd == "-h" || a.cmd == "--help" || a.cmd == "help" || a.has("help") || a.has("h")) {
    std::fputs(kUsage, stdout);
    return 0;
  }
  if (a.cmd == "spectrogram" || a.cmd == "measure" || a.cmd == "locate") {
    bool on = gpu::init(!a.has("no-gpu"), (int)a.num("gpu-device", 0));
    if (gpu::compiled()) log("gpu: %s", on ? gpu::info().c_str() : gpu::info().c_str());
  }
  if (a.cmd == "spectrogram") return cmd_spectrogram(a);
  if (a.cmd == "measure") return cmd_measure(a);
  if (a.cmd == "locate") return cmd_locate(a);
  if (a.cmd == "ingest") return cmd_ingest(a);
  if (a.cmd == "simulate") return cmd_simulate(a);
  if (a.cmd == "update") return cmd_update(a);
  if (a.cmd == "mcp") return cmd_mcp(a);
  die("unknown command '%s' (try --help)", a.cmd.c_str());
}
