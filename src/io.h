// Small self-contained I/O helpers: PNG writer, minimal JSON DOM, gzip-aware
// file reading.
#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace pl {

struct RGB {
  uint8_t r, g, b;
};

struct Image {
  int w = 0, h = 0;
  std::vector<RGB> px;
  Image(int W, int H) : w(W), h(H), px((size_t)W * H, RGB{0, 0, 0}) {}
  void set(int x, int y, RGB c) {
    if (x >= 0 && y >= 0 && x < w && y < h) px[(size_t)y * w + x] = c;
  }
  void line(int x0, int y0, int x1, int y1, RGB c);
  void dot(int x, int y, int r, RGB c);
};

void write_png(const std::string& path, const Image& img);
RGB colormap(float t);  // t in [0,1], perceptual "inferno"-like ramp

// Reads a whole file; transparently gunzips if the content is gzip-compressed.
bool read_file_maybe_gz(const std::string& path, std::string& out);

struct JMember;

struct JVal {
  enum Type : uint8_t { Null, Bool, Num, Str, Arr, Obj } t = Null;
  bool b = false;
  double n = 0;
  std::string s;
  std::vector<JVal> a;
  // Members in document order. A vector of std::pair<std::string, JVal> would
  // instantiate the pair while JVal is incomplete, which Clang rejects.
  std::vector<JMember> o;
  const JVal* get(const char* key) const;
  bool is_num() const { return t == Num; }
};
struct JMember {
  std::string first;
  JVal second;
};

bool parse_json(const std::string& text, JVal& out, std::string* err = nullptr);
std::string json_escape(const std::string& s);

}  // namespace pl
