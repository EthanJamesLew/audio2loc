// Outputs: GeoJSON layers, a ready-to-open QGIS project, spectrogram PNGs with model overlays.
#pragma once
#include <string>
#include <utility>
#include <vector>

#include "dsp.h"
#include "io.h"

namespace pl {

struct PointFeature {
  double lat, lon;
  std::string props;  // JSON members without braces, e.g. "\"score\":1.2,\"id\":\"x\""
};
struct LineFeature {
  std::vector<std::pair<double, double>> lonlat;
  std::string props;
};
void write_points_geojson(const std::string& path, const std::vector<PointFeature>& f);
void write_lines_geojson(const std::string& path, const std::vector<LineFeature>& f);

// Incremental point-layer writer (constant memory for arbitrarily large layers).
class GeoJsonWriter {
 public:
  ~GeoJsonWriter() { close(); }
  void open(const std::string& path);
  void add(double lat, double lon, const std::string& props);
  void close();

 private:
  FILE* f_ = nullptr;
  size_t n_ = 0;
};

// Layer styles understood by both the XML writer and the native QGIS exporter.
struct LayerSpec {
  std::string name, file;  // file relative to the project directory
  enum Kind { Points, Lines } kind = Points;
  enum Style { Single, Graduated, Heatmap } style = Single;
  std::string field;                // graduated / heatmap weight field
  std::vector<double> breaks;       // graduated class breaks (ascending, n+1 values)
  std::string color = "230,30,30";  // single-symbol color "r,g,b"
  double size = 2.0;                // marker size / line width (mm)
  std::string marker = "circle";
  bool visible = true;
};

struct ProjectSpec {
  std::string title;
  std::vector<LayerSpec> layers;  // drawn top to bottom
  double minlon, minlat, maxlon, maxlat;
  bool osm_basemap = true;
};

void write_qgs_project(const std::string& path, const ProjectSpec& p);
// Native QGIS C++ API export (.qgz). Returns false if this build has no QGIS support.
bool write_qgis_native(const std::string& path, const std::string& dir, const ProjectSpec& p);
bool qgis_native_available();

struct Curve {
  std::vector<std::pair<double, double>> tf;  // (audio time, Hz)
  RGB color;
};
Image render_spectrogram(const Spectrogram& white, double t0, double t1, double f0, double f1,
                         const std::vector<Curve>& curves, int yscale = 1);

std::vector<double> quantile_breaks(std::vector<double> v, int n);

}  // namespace pl
