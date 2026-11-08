// Model Context Protocol server (stdio, JSON-RPC 2.0, newline-delimited)
// exposing audio2loc to AI agents. Each tool call runs the audio2loc binary
// itself as a subprocess so long searches or crashes never take the server
// down; images (spectrograms, model overlays) are returned as MCP image content
// so a multimodal model can find and box Doppler curves itself.
#include "mcp.h"
#include "version.h"

#include <unistd.h>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "io.h"
#include "update.h"
#include "util.h"

namespace fs = std::filesystem;

namespace pl {
namespace {

std::string g_self, g_work, g_db;

std::string q(const std::string& s) {
  return "\"" + json_escape(s) + "\"";
}

std::string shq(const std::string& s) {
  std::string o = "'";
  for (char c : s)
    o += c == '\'' ? std::string("'\\''") : std::string(1, c);
  return o + "'";
}

std::string b64(const std::string& in) {
  static const char* t =
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  std::string o;
  o.reserve((in.size() + 2) / 3 * 4);
  size_t i = 0;
  for (; i + 2 < in.size(); i += 3) {
    uint32_t v = ((unsigned char)in[i] << 16) |
                 ((unsigned char)in[i + 1] << 8) | (unsigned char)in[i + 2];
    o += t[v >> 18], o += t[(v >> 12) & 63], o += t[(v >> 6) & 63],
        o += t[v & 63];
  }
  if (i < in.size()) {
    uint32_t v = (unsigned char)in[i] << 16;
    if (i + 1 < in.size()) v |= (unsigned char)in[i + 1] << 8;
    o += t[v >> 18], o += t[(v >> 12) & 63];
    o += i + 1 < in.size() ? t[(v >> 6) & 63] : '=';
    o += '=';
  }
  return o;
}

std::string id_json(const JVal* id) {
  if (!id) return "null";
  if (id->t == JVal::Str) return q(id->s);
  if (id->t == JVal::Num) {
    char b[64];
    std::snprintf(b, sizeof b, "%.17g", id->n);
    return b;
  }
  return "null";
}

void send(const std::string& msg) {
  std::cout << msg << "\n" << std::flush;
}

// ---------------------------------------------------------------- tools ---
struct Content {
  std::vector<std::string> items;  // serialized content blocks
  bool error = false;
  void text(const std::string& s) {
    items.push_back("{\"type\":\"text\",\"text\":" + q(s) + "}");
  }
  void png(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return;
    std::string d((std::istreambuf_iterator<char>(f)),
                  std::istreambuf_iterator<char>());
    items.push_back(
        "{\"type\":\"image\",\"mimeType\":\"image/png\",\"data\":\"" + b64(d) +
        "\"}");
  }
};

int run(const std::vector<std::string>& args, std::string& out) {
  std::string cmd = shq(g_self);
  for (auto& a : args)
    cmd += " " + shq(a);
  cmd += " 2>&1 </dev/null";
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) return -1;
  char buf[8192];
  size_t n;
  std::string raw;
  while ((n = std::fread(buf, 1, sizeof buf, p)) > 0)
    raw.append(buf, n);
  int rc = pclose(p);
  // drop carriage-return progress lines
  out.clear();
  std::string line;
  std::istringstream is(raw);
  while (std::getline(is, line)) {
    size_t cr = line.rfind('\r');
    if (cr != std::string::npos) line = line.substr(cr + 1);
    out += line + "\n";
  }
  if (out.size() > 24000)
    out = out.substr(0, 4000) + "\n...[truncated]...\n" +
          out.substr(out.size() - 18000);
  return WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
}

std::string num_or(const JVal* v, const char* def = nullptr) {
  if (v && v->t == JVal::Num) {
    char b[64];
    std::snprintf(b, sizeof b, "%.10g", v->n);
    return b;
  }
  if (v && v->t == JVal::Str) return v->s;
  return def ? def : "";
}

// Local path, or download http(s) URLs (yt-dlp for pages, curl for direct
// media) into the work dir.
bool resolve_audio(const JVal& args, std::string& path, Content& c) {
  const JVal* a = args.get("audio");
  if (!a || a->t != JVal::Str || a->s.empty()) {
    c.text(
        "missing 'audio' (path inside the server's filesystem, or an http(s) "
        "URL)");
    return false;
  }
  if (a->s.rfind("http://", 0) != 0 && a->s.rfind("https://", 0) != 0) {
    path = a->s;
    if (!fs::exists(path) && fs::exists(g_work + "/" + path))
      path = g_work + "/" + path;
    if (!fs::exists(path)) {
      c.text("audio file not found: " + a->s);
      return false;
    }
    return true;
  }
  std::string dir = g_work + "/media";
  fs::create_directories(dir);
  std::string h = std::to_string(std::hash<std::string>{}(a->s));
  for (auto& e : fs::directory_iterator(dir))
    if (e.path().filename().string().rfind(h + ".", 0) == 0) {
      path = e.path().string();
      return true;
    }
  // stdout is the JSON-RPC channel: everything from the downloaders goes to
  // stderr.
  std::string cmd = "{ yt-dlp -q --no-playlist -f 'bestaudio/best' -o " +
                    shq(dir + "/" + h + ".%(ext)s") + " " + shq(a->s) +
                    " || curl -fsSL -o " + shq(dir + "/" + h + ".media") + " " +
                    shq(a->s) + "; } 1>&2 </dev/null";
  if (std::system(cmd.c_str()) != 0) {
    c.text("could not download " + a->s + " (yt-dlp and curl both failed)");
    return false;
  }
  for (auto& e : fs::directory_iterator(dir))
    if (e.path().filename().string().rfind(h + ".", 0) == 0)
      path = e.path().string();
  return !path.empty();
}

bool event_args(const JVal& args, std::vector<std::string>& out, Content& c) {
  const JVal* ev = args.get("events");
  if (!ev || ev->t != JVal::Arr || ev->a.empty()) {
    c.text(
        "missing 'events': array of {t0,t1,fmin,fmax,...} boxes around Doppler "
        "curves");
    return false;
  }
  for (auto& e : ev->a) {
    std::string s;
    for (const char* k : {"t0", "t1", "fmin", "fmax", "tc", "tc_tol"}) {
      std::string v = num_or(e.get(k));
      if (!v.empty()) s += std::string(s.empty() ? "" : ",") + k + "=" + v;
    }
    auto arr = [&](const char* k, char sep) {
      const JVal* v = e.get(k);
      if (!v || v->t != JVal::Arr || v->a.empty()) return;
      s += std::string(",") + k + "=";
      for (size_t i = 0; i < v->a.size(); ++i)
        s += (i ? std::string(1, sep) : "") + num_or(&v->a[i]);
    };
    arr("mult", ':');
    arr("v", ':');
    arr("d", ':');
    arr("f0", ':');
    if (auto* v = e.get("side"); v && v->t == JVal::Str) s += ",side=" + v->s;
    if (auto* v = e.get("score"); v && v->t == JVal::Str) s += ",score=" + v->s;
    out.push_back("--event");
    out.push_back(s);
  }
  return true;
}

std::string new_run_dir(const char* kind) {
  static std::atomic<int> seq{0};
  std::string d = g_work + "/runs/" + utc_date((double)std::time(nullptr)) +
                  "_" + kind + "_" + std::to_string((long)std::time(nullptr)) +
                  "_" + std::to_string(seq++);
  fs::create_directories(d);
  return d;
}

void tool_spectrogram(const JVal& a, Content& c) {
  std::string audio;
  if (!resolve_audio(a, audio, c)) return void(c.error = true);
  std::string dir = new_run_dir("spec"), png = dir + "/spectrogram.png", out;
  std::vector<std::string> args = {"spectrogram", "--audio", audio, "-o", png};
  for (const char* k : {"t0", "t1", "fmin", "fmax", "hop", "win"}) {
    std::string v = num_or(a.get(k));
    if (!v.empty()) args.push_back(std::string("--") + k), args.push_back(v);
  }
  if (!a.get("fmin")) args.push_back("--fmin"), args.push_back("40");
  if (!a.get("fmax")) args.push_back("--fmax"), args.push_back("400");
  int rc = run(args, out);
  c.text(out + (rc ? ""
                   : "Image: x = time (left->right), y = frequency (high at "
                     "top). Doppler curves of propeller aircraft "
                     "look like falling S-shaped tones (often with parallel "
                     "harmonics) lasting ~1-3 min. Use the scale "
                     "line above to convert pixels to seconds/Hz, then box ONE "
                     "tone of each curve for measure_doppler."));
  c.error = rc != 0;
  if (!rc) c.png(png);
}

void tool_measure(const JVal& a, Content& c) {
  std::string audio;
  std::vector<std::string> args = {"measure"};
  if (!resolve_audio(a, audio, c) || !event_args(a, args, c))
    return void(c.error = true);
  std::string dir = new_run_dir("measure"), out;
  args.insert(args.end(), {"--audio", audio, "--out", dir});
  int rc = run(args, out);
  c.text(out);
  c.error = rc != 0;
  for (int k = 1; k <= 8; ++k)
    c.png(dir + "/measure_event" + std::to_string(k) + ".png");
}

void tool_locate(const JVal& a, Content& c) {
  std::string audio;
  std::vector<std::string> args = {"locate"};
  if (!resolve_audio(a, audio, c) || !event_args(a, args, c))
    return void(c.error = true);
  std::string dir = new_run_dir("locate"), out;
  args.insert(args.end(), {"--audio", audio, "--out", dir, "--overlays", "2"});
  bool own_adsb = false;
  if (auto* x = a.get("extra_args"); x && x->t == JVal::Arr)
    for (auto& s : x->a)
      own_adsb |= s.t == JVal::Str && s.s == "--adsb";
  if (!own_adsb) args.insert(args.end(), {"--db", g_db, "--auto-update"});
  bool has_time = false;
  for (const char* k : {"start", "end", "rec_start", "rec_start_tol"}) {
    std::string v = num_or(a.get(k));
    if (v.empty()) continue;
    std::string flag = k;
    for (auto& ch : flag)
      if (ch == '_') ch = '-';
    args.push_back("--" + flag), args.push_back(v);
    has_time |= std::string(k) == "start" || std::string(k) == "rec_start";
  }
  if (!has_time) {
    c.text(
        "locate needs a time window: 'start' (+ optional 'end') or 'rec_start' "
        "(ISO 8601 UTC)");
    return void(c.error = true);
  }
  if (auto* b = a.get("bbox"); b && b->t == JVal::Arr && b->a.size() == 4) {
    std::string s;
    for (int i = 0; i < 4; ++i)
      s += (i ? "," : "") + num_or(&b->a[i]);
    args.push_back("--bbox"), args.push_back(s);
  }
  if (auto* v = a.get("include_jets"); v && v->t == JVal::Bool && v->b)
    args.push_back("--include-jets");
  if (auto* x = a.get("extra_args"); x && x->t == JVal::Arr)
    for (auto& s : x->a)
      if (s.t == JVal::Str) args.push_back(s.s);
  int rc = run(args, out);
  std::string res;
  read_file_maybe_gz(dir + "/result.json", res);
  c.text(out + (res.empty() ? "" : "\nresult.json:\n" + res) +
         "\nOutput directory (GeoJSON layers + QGIS project audio2loc.qgs): " +
         dir);
  c.error = rc != 0;
  for (int k = 1; k <= 8; ++k)
    c.png(dir + "/overlay_rank1_event" + std::to_string(k) + ".png");
}

void tool_update(const JVal& a, Content& c) {
  std::vector<std::string> args = {"update", "--db", g_db};
  if (auto* d = a.get("dates"); d && d->t == JVal::Arr)
    for (auto& x : d->a)
      if (x.t == JVal::Str) args.push_back("--date"), args.push_back(x.s);
  if (auto* v = a.get("latest"); v && v->t == JVal::Bool && v->b)
    args.push_back("--latest");
  if (args.size() == 3) args.push_back("--latest");
  if (auto* b = a.get("bbox"); b && b->t == JVal::Arr && b->a.size() == 4) {
    std::string s;
    for (int i = 0; i < 4; ++i)
      s += (i ? "," : "") + num_or(&b->a[i]);
    args.push_back("--bbox"), args.push_back(s);
  }
  std::string out;
  int rc = run(args, out);
  c.text(out);
  c.error = rc != 0;
}

void tool_status(const JVal&, Content& c) {
  std::string out;
  run({"update", "--db", g_db, "--list"}, out);
  c.text(out +
         "Newest adsb.lol day is normally yesterday (UTC); dumps publish "
         "~03:30 UTC the next day.");
}

const char* kTools =
    R"JSON([{"name":"spectrogram","description":"Render a whitened spectrogram PNG of an audio/video file (or URL) to find aircraft Doppler curves. Start broad (whole clip, 40-400 Hz, hop 0.5 for long clips), then zoom (shorter t0..t1, narrower band) to read box coordinates precisely.","inputSchema":{"type":"object","properties":{"audio":{"type":"string","description":"Path in the server filesystem (e.g. /work/clip.mp4) or an http(s) URL (downloaded with yt-dlp)"},"t0":{"type":"number","description":"start (s)"},"t1":{"type":"number","description":"end (s)"},"fmin":{"type":"number","description":"Hz, default 40"},"fmax":{"type":"number","description":"Hz, default 400"},"hop":{"type":"number","description":"seconds per column, default 0.1"}},"required":["audio"]}},{"name":"measure_doppler","description":"Fit boxed Doppler curves: returns aircraft speed, closest-approach slant distance, CPA time and source frequency (with uncertainty ranges) plus overlay PNGs (green: traced ridge, cyan: fit). Check the overlays before running locate.","inputSchema":{"type":"object","properties":{"audio":{"type":"string","description":"Path in the server filesystem (e.g. /work/clip.mp4) or an http(s) URL (downloaded with yt-dlp)"},"events":{"type":"array","description":"one box per Doppler curve, around ONE tone of it","items":{"type":"object","properties":{"t0":{"type":"number","description":"box start (s)"},"t1":{"type":"number","description":"box end (s)"},"fmin":{"type":"number","description":"box low freq (Hz)"},"fmax":{"type":"number","description":"box high freq (Hz)"},"mult":{"type":"array","items":{"type":"number"},"description":"harmonics of the boxed tone used in scoring, e.g. [1,2,3]"},"side":{"type":"string","enum":["both","before","after"],"description":"which side of CPA to fit (use if one side is distorted)"},"score":{"type":"string","enum":["brightness","ridge"]},"v":{"type":"array","items":{"type":"number"},"description":"override speed range [min,max] m/s"},"d":{"type":"array","items":{"type":"number"},"description":"override slant distance range [min,max] m"},"f0":{"type":"array","items":{"type":"number"},"description":"override source frequency range Hz"},"tc":{"type":"number","description":"override CPA time (s)"},"tc_tol":{"type":"number"}},"required":["t0","t1","fmin","fmax"]}}},"required":["audio","events"]}},{"name":"locate","description":"Full geolocation search: matches the Doppler events against ADS-B history (auto-downloads needed days), applies time-offset and silence filters, models every candidate's exact Doppler curve from real trajectories, and ranks locations. Returns ranked coordinates, flights involved, implied recording start time, overlay PNGs of the winning model and an output dir with GeoJSON + QGIS project. Needs 2+ events, or 1 event plus a tight rec_start, for a unique fix. Can take many minutes for large areas.","inputSchema":{"type":"object","properties":{"audio":{"type":"string","description":"Path in the server filesystem (e.g. /work/clip.mp4) or an http(s) URL (downloaded with yt-dlp)"},"events":{"type":"array","description":"one box per Doppler curve, around ONE tone of it","items":{"type":"object","properties":{"t0":{"type":"number","description":"box start (s)"},"t1":{"type":"number","description":"box end (s)"},"fmin":{"type":"number","description":"box low freq (Hz)"},"fmax":{"type":"number","description":"box high freq (Hz)"},"mult":{"type":"array","items":{"type":"number"},"description":"harmonics of the boxed tone used in scoring, e.g. [1,2,3]"},"side":{"type":"string","enum":["both","before","after"],"description":"which side of CPA to fit (use if one side is distorted)"},"score":{"type":"string","enum":["brightness","ridge"]},"v":{"type":"array","items":{"type":"number"},"description":"override speed range [min,max] m/s"},"d":{"type":"array","items":{"type":"number"},"description":"override slant distance range [min,max] m"},"f0":{"type":"array","items":{"type":"number"},"description":"override source frequency range Hz"},"tc":{"type":"number","description":"override CPA time (s)"},"tc_tol":{"type":"number"}},"required":["t0","t1","fmin","fmax"]}},"start":{"type":"string","description":"ISO 8601 UTC start of the window in which the recording was made"},"end":{"type":"string","description":"ISO 8601 UTC end (default start + 24 h)"},"rec_start":{"type":"string","description":"approximate absolute UTC time of audio t=0, if known"},"rec_start_tol":{"type":"number","description":"uncertainty of rec_start in s (default 300)"},"bbox":{"type":"array","items":{"type":"number"},"description":"[W,S,E,N] degrees; strongly recommended for speed"},"include_jets":{"type":"boolean"},"extra_args":{"type":"array","items":{"type":"string"},"description":"additional audio2loc locate CLI flags"}},"required":["audio","events"]}},{"name":"update_database","description":"Download and index daily global ADS-B history from adsb.lol into the local database (one UTC day is a ~3-4 GB streamed download, ~10-20 min). The newest available day is normally yesterday (UTC).","inputSchema":{"type":"object","properties":{"dates":{"type":"array","items":{"type":"string"},"description":"YYYY-MM-DD (UTC)"},"latest":{"type":"boolean","description":"fetch the newest published day"},"bbox":{"type":"array","items":{"type":"number"},"description":"[W,S,E,N]: keep only this region (smaller, faster to search)"}}}},{"name":"database_status","description":"List the UTC days available in the local ADS-B database.","inputSchema":{"type":"object","properties":{}}}])JSON";

const char* kInstructions =
    "audio2loc geolocates a recording from the Doppler curves of propeller "
    "aircraft flying overhead, matched "
    "against ADS-B history. Workflow: (1) spectrogram over the whole clip "
    "40-400 Hz; look for falling S-shaped "
    "tones lasting 1-3 minutes, often inaudible. (2) zoom into each and box "
    "ONE clear tone (t0,t1,fmin,fmax); "
    "add mult [1,2,3] if harmonics are visible. (3) measure_doppler and "
    "inspect the overlays: the cyan fit should "
    "follow the green ridge. (4) locate with the date window (and bbox / "
    "rec_start when known). Treat a fix as "
    "credible only when the overlays show the modelled curve covering the "
    "observed curve and #1 clearly beats #2. "
    "Intended for verification/OSINT of material whose subject consents or "
    "that is of public interest; do not use "
    "it to find private individuals.";

}  // namespace

int run_mcp_server(const std::string& self, const std::string& work,
                   const std::string& db) {
  g_self = self;
  g_work = work;
  g_db = db;
  fs::create_directories(g_work);
  std::fprintf(stderr, "audio2loc MCP server (stdio). work=%s db=%s\n",
               g_work.c_str(), g_db.c_str());
  std::string line;
  while (std::getline(std::cin, line)) {
    if (trim(line).empty()) continue;
    JVal m;
    if (!parse_json(line, m) || m.t != JVal::Obj) {
      send(
          "{\"jsonrpc\":\"2.0\",\"id\":null,\"error\":{\"code\":-32700,"
          "\"message\":\"parse error\"}}");
      continue;
    }
    const JVal* method = m.get("method");
    const JVal* id = m.get("id");
    if (!method || method->t != JVal::Str)
      continue;  // a response to us; ignore
    std::string meth = method->s;
    if (!id) continue;  // notifications (initialized, cancelled, ...)
    std::string rid = id_json(id);
    const JVal* params = m.get("params");
    JVal empty;
    empty.t = JVal::Obj;
    if (!params) params = &empty;
    if (meth == "initialize") {
      std::string pv = "2025-06-18";
      if (auto* v = params->get("protocolVersion");
          v && v->t == JVal::Str &&
          (v->s == "2024-11-05" || v->s == "2025-03-26" ||
           v->s == "2025-06-18"))
        pv = v->s;
      send("{\"jsonrpc\":\"2.0\",\"id\":" + rid +
           ",\"result\":{\"protocolVersion\":" + q(pv) +
           ",\"capabilities\":{\"tools\":{\"listChanged\":false}},"
           "\"serverInfo\":{\"name\":\"audio2loc\",\"version\":" +
           q(version_string()) + "},\"instructions\":" + q(kInstructions) +
           "}}");
    } else if (meth == "ping") {
      send("{\"jsonrpc\":\"2.0\",\"id\":" + rid + ",\"result\":{}}");
    } else if (meth == "tools/list") {
      send("{\"jsonrpc\":\"2.0\",\"id\":" + rid +
           ",\"result\":{\"tools\":" + kTools + "}}");
    } else if (meth == "tools/call") {
      const JVal* name = params->get("name");
      const JVal* args = params->get("arguments");
      if (!args) args = &empty;
      Content c;
      std::string n = name && name->t == JVal::Str ? name->s : "";
      if (n == "spectrogram")
        tool_spectrogram(*args, c);
      else if (n == "measure_doppler")
        tool_measure(*args, c);
      else if (n == "locate")
        tool_locate(*args, c);
      else if (n == "update_database")
        tool_update(*args, c);
      else if (n == "database_status")
        tool_status(*args, c);
      else {
        send("{\"jsonrpc\":\"2.0\",\"id\":" + rid +
             ",\"error\":{\"code\":-32602,\"message\":" +
             q("unknown tool " + n) + "}}");
        continue;
      }
      std::string body;
      for (size_t i = 0; i < c.items.size(); ++i)
        body += (i ? "," : "") + c.items[i];
      send("{\"jsonrpc\":\"2.0\",\"id\":" + rid + ",\"result\":{\"content\":[" +
           body + "],\"isError\":" + (c.error ? "true" : "false") + "}}");
    } else {
      send("{\"jsonrpc\":\"2.0\",\"id\":" + rid +
           ",\"error\":{\"code\":-32601,\"message\":" +
           q("method not found: " + meth) + "}}");
    }
  }
  return 0;
}

}  // namespace pl
