#include "io.h"

#include <zlib.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "util.h"

namespace pl {

// ---------------------------------------------------------------- PNG ---
void Image::line(int x0, int y0, int x1, int y1, RGB c) {
  int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx + dy;
  for (int guard = 0; guard < 100000; ++guard) {
    set(x0, y0, c);
    if (x0 == x1 && y0 == y1) break;
    int e2 = 2 * err;
    if (e2 >= dy) err += dy, x0 += sx;
    if (e2 <= dx) err += dx, y0 += sy;
  }
}

void Image::dot(int x, int y, int r, RGB c) {
  for (int dy = -r; dy <= r; ++dy)
    for (int dx = -r; dx <= r; ++dx)
      if (dx * dx + dy * dy <= r * r) set(x + dx, y + dy, c);
}

static void put32(std::vector<unsigned char>& v, uint32_t x) {
  v.push_back(x >> 24);
  v.push_back(x >> 16);
  v.push_back(x >> 8);
  v.push_back(x);
}

static void chunk(FILE* f, const char* type, const std::vector<unsigned char>& data) {
  std::vector<unsigned char> buf;
  put32(buf, (uint32_t)data.size());
  buf.insert(buf.end(), type, type + 4);
  buf.insert(buf.end(), data.begin(), data.end());
  uLong crc = crc32(0, buf.data() + 4, (uInt)(buf.size() - 4));
  put32(buf, (uint32_t)crc);
  std::fwrite(buf.data(), 1, buf.size(), f);
}

void write_png(const std::string& path, const Image& img) {
  FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) die("cannot write %s", path.c_str());
  static const unsigned char sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
  std::fwrite(sig, 1, 8, f);
  std::vector<unsigned char> ihdr;
  put32(ihdr, img.w);
  put32(ihdr, img.h);
  ihdr.insert(ihdr.end(), {8, 2, 0, 0, 0});
  chunk(f, "IHDR", ihdr);
  std::vector<unsigned char> raw;
  raw.reserve((size_t)img.h * (img.w * 3 + 1));
  for (int y = 0; y < img.h; ++y) {
    raw.push_back(0);
    for (int x = 0; x < img.w; ++x) {
      const RGB& p = img.px[(size_t)y * img.w + x];
      raw.insert(raw.end(), {p.r, p.g, p.b});
    }
  }
  uLongf clen = compressBound((uLong)raw.size());
  std::vector<unsigned char> comp(clen);
  compress2(comp.data(), &clen, raw.data(), (uLong)raw.size(), 6);
  comp.resize(clen);
  chunk(f, "IDAT", comp);
  chunk(f, "IEND", {});
  std::fclose(f);
}

RGB colormap(float t) {
  static const float k[][3] = {{0, 0, 4},      {40, 11, 84},   {101, 21, 110}, {159, 42, 99},
                               {212, 72, 66},  {245, 125, 21}, {250, 193, 39}, {252, 255, 164}};
  t = std::clamp(t, 0.f, 1.f) * 7.f;
  int i = std::min(6, (int)t);
  float u = t - i;
  return RGB{(uint8_t)(k[i][0] + (k[i + 1][0] - k[i][0]) * u), (uint8_t)(k[i][1] + (k[i + 1][1] - k[i][1]) * u),
             (uint8_t)(k[i][2] + (k[i + 1][2] - k[i][2]) * u)};
}

bool read_file_maybe_gz(const std::string& path, std::string& out) {
  gzFile g = gzopen(path.c_str(), "rb");
  if (!g) return false;
  gzbuffer(g, 1 << 17);
  out.clear();
  char buf[1 << 16];
  int n;
  while ((n = gzread(g, buf, sizeof buf)) > 0) out.append(buf, n);
  gzclose(g);
  return n >= 0;
}

// --------------------------------------------------------------- JSON ---
const JVal* JVal::get(const char* key) const {
  if (t != Obj) return nullptr;
  for (auto& kv : o)
    if (kv.first == key) return &kv.second;
  return nullptr;
}

namespace {
struct JP {
  const char* p;
  const char* e;
  std::string err;
  void ws() {
    while (p < e && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')) ++p;
  }
  bool str(std::string& s) {
    ++p;  // opening quote
    while (p < e && *p != '"') {
      if (*p == '\\' && p + 1 < e) {
        ++p;
        switch (*p) {
          case 'n': s += '\n'; break;
          case 't': s += '\t'; break;
          case 'r': s += '\r'; break;
          case 'b': s += '\b'; break;
          case 'f': s += '\f'; break;
          case 'u': {
            if (e - p < 5) return false;
            unsigned cp = std::strtoul(std::string(p + 1, p + 5).c_str(), nullptr, 16);
            if (cp < 0x80) s += (char)cp;
            else if (cp < 0x800) s += (char)(0xC0 | (cp >> 6)), s += (char)(0x80 | (cp & 0x3F));
            else s += (char)(0xE0 | (cp >> 12)), s += (char)(0x80 | ((cp >> 6) & 0x3F)), s += (char)(0x80 | (cp & 0x3F));
            p += 4;
            break;
          }
          default: s += *p;
        }
        ++p;
      } else {
        s += *p++;
      }
    }
    if (p >= e) return false;
    ++p;
    return true;
  }
  bool val(JVal& v, int depth) {
    if (depth > 200) return false;
    ws();
    if (p >= e) return false;
    char c = *p;
    if (c == '{') {
      v.t = JVal::Obj;
      ++p;
      ws();
      if (p < e && *p == '}') return ++p, true;
      for (;;) {
        ws();
        if (p >= e || *p != '"') return false;
        std::string k;
        if (!str(k)) return false;
        ws();
        if (p >= e || *p != ':') return false;
        ++p;
        v.o.emplace_back(std::move(k), JVal{});
        if (!val(v.o.back().second, depth + 1)) return false;
        ws();
        if (p < e && *p == ',') { ++p; continue; }
        if (p < e && *p == '}') return ++p, true;
        return false;
      }
    }
    if (c == '[') {
      v.t = JVal::Arr;
      ++p;
      ws();
      if (p < e && *p == ']') return ++p, true;
      for (;;) {
        v.a.emplace_back();
        if (!val(v.a.back(), depth + 1)) return false;
        ws();
        if (p < e && *p == ',') { ++p; continue; }
        if (p < e && *p == ']') return ++p, true;
        return false;
      }
    }
    if (c == '"') {
      v.t = JVal::Str;
      return str(v.s);
    }
    if (c == 't' && e - p >= 4 && !std::memcmp(p, "true", 4)) return v.t = JVal::Bool, v.b = true, p += 4, true;
    if (c == 'f' && e - p >= 5 && !std::memcmp(p, "false", 5)) return v.t = JVal::Bool, p += 5, true;
    if (c == 'n' && e - p >= 4 && !std::memcmp(p, "null", 4)) return v.t = JVal::Null, p += 4, true;
    char* end;
    v.n = std::strtod(p, &end);
    if (end == p) return false;
    v.t = JVal::Num;
    p = end;
    return true;
  }
};
}  // namespace

bool parse_json(const std::string& text, JVal& out, std::string* err) {
  JP jp{text.data(), text.data() + text.size(), {}};
  bool ok = jp.val(out, 0);
  if (!ok && err) *err = "parse error at offset " + std::to_string(jp.p - text.data());
  return ok;
}

std::string json_escape(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\', o += c;
    else if ((unsigned char)c < 0x20) {
      char b[8];
      std::snprintf(b, sizeof b, "\\u%04x", c);
      o += b;
    } else o += c;
  }
  return o;
}

}  // namespace pl
