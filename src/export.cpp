#include "export.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "geo.h"
#include "util.h"

namespace pl {

void write_points_geojson(const std::string& path,
                          const std::vector<PointFeature>& f) {
  FILE* o = std::fopen(path.c_str(), "w");
  if (!o) die("cannot write %s", path.c_str());
  std::fprintf(o,
               "{\"type\":\"FeatureCollection\",\"crs\":{\"type\":\"name\","
               "\"properties\":{\"name\":\"urn:ogc:def:crs:OGC:1.3:CRS84\"}},"
               "\"features\":[\n");
  for (size_t i = 0; i < f.size(); ++i)
    std::fprintf(o,
                 "{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\","
                 "\"coordinates\":[%.7f,%.7f]},\"properties\":{%s}}%s\n",
                 f[i].lon, f[i].lat, f[i].props.c_str(),
                 i + 1 < f.size() ? "," : "");
  std::fprintf(o, "]}\n");
  std::fclose(o);
}

void write_lines_geojson(const std::string& path,
                         const std::vector<LineFeature>& f) {
  FILE* o = std::fopen(path.c_str(), "w");
  if (!o) die("cannot write %s", path.c_str());
  std::fprintf(o, "{\"type\":\"FeatureCollection\",\"features\":[\n");
  for (size_t i = 0; i < f.size(); ++i) {
    std::fprintf(o,
                 "{\"type\":\"Feature\",\"geometry\":{\"type\":\"LineString\","
                 "\"coordinates\":[");
    for (size_t k = 0; k < f[i].lonlat.size(); ++k)
      std::fprintf(o, "%s[%.7f,%.7f]", k ? "," : "", f[i].lonlat[k].first,
                   f[i].lonlat[k].second);
    std::fprintf(o, "]},\"properties\":{%s}}%s\n", f[i].props.c_str(),
                 i + 1 < f.size() ? "," : "");
  }
  std::fprintf(o, "]}\n");
  std::fclose(o);
}

void GeoJsonWriter::open(const std::string& path) {
  close();
  f_ = std::fopen(path.c_str(), "w");
  if (!f_) die("cannot write %s", path.c_str());
  n_ = 0;
  std::fprintf(f_, "{\"type\":\"FeatureCollection\",\"features\":[\n");
}

void GeoJsonWriter::add(double lat, double lon, const std::string& props) {
  if (!f_) return;
  std::fprintf(f_,
               "%s{\"type\":\"Feature\",\"geometry\":{\"type\":\"Point\","
               "\"coordinates\":[%.7f,%.7f]},\"properties\":{%s}}",
               n_++ ? ",\n" : "", lon, lat, props.c_str());
}

void GeoJsonWriter::close() {
  if (!f_) return;
  std::fprintf(f_, "\n]}\n");
  std::fclose(f_);
  f_ = nullptr;
}

std::vector<double> quantile_breaks(std::vector<double> v, int n) {
  std::vector<double> b;
  if (v.empty()) return {0, 1};
  std::sort(v.begin(), v.end());
  for (int i = 0; i <= n; ++i) {
    double x = v[std::min(v.size() - 1,
                          (size_t)std::lround((double)i / n * (v.size() - 1)))];
    if (b.empty() || x > b.back() + 1e-9) b.push_back(x);
  }
  if (b.size() < 2) b.push_back(b.back() + 1e-6);
  b.front() -= 1e-6;
  b.back() += 1e-6;
  return b;
}

// --------------------------------------------------------------- QGIS ---
static std::string xml(const std::string& s) {
  std::string o;
  for (char c : s) {
    switch (c) {
      case '&':
        o += "&amp;";
        break;
      case '<':
        o += "&lt;";
        break;
      case '>':
        o += "&gt;";
        break;
      case '"':
        o += "&quot;";
        break;
      default:
        o += c;
    }
  }
  return o;
}

static const char* kSrs4326 =
    "<spatialrefsys nativeFormat=\"Wkt\"><wkt>GEOGCRS[\"WGS 84\",DATUM[\"World "
    "Geodetic System 1984\",ELLIPSOID[\"WGS "
    "84\",6378137,298.257223563,LENGTHUNIT[\"metre\",1]]],PRIMEM[\"Greenwich\","
    "0,ANGLEUNIT[\"degree\",0.0174532925199433]],"
    "CS[ellipsoidal,2],AXIS[\"geodetic latitude "
    "(Lat)\",north,ORDER[1],ANGLEUNIT[\"degree\",0.0174532925199433]],"
    "AXIS[\"geodetic longitude "
    "(Lon)\",east,ORDER[2],ANGLEUNIT[\"degree\",0.0174532925199433]],ID["
    "\"EPSG\",4326]]</wkt>"
    "<proj4>+proj=longlat +datum=WGS84 "
    "+no_defs</proj4><srsid>3452</srsid><srid>4326</srid><authid>EPSG:4326</"
    "authid>"
    "<description>WGS "
    "84</description><projectionacronym>longlat</"
    "projectionacronym><ellipsoidacronym>EPSG:7030</ellipsoidacronym>"
    "<geographicflag>true</geographicflag></spatialrefsys>";
static const char* kSrs3857 =
    "<spatialrefsys nativeFormat=\"Wkt\"><proj4>+proj=merc +a=6378137 "
    "+b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 "
    "+units=m +nadgrids=@null +wktext "
    "+no_defs</proj4><srsid>3857</srsid><srid>3857</srid><authid>EPSG:3857</"
    "authid>"
    "<description>WGS 84 / "
    "Pseudo-Mercator</description><projectionacronym>merc</projectionacronym>"
    "<ellipsoidacronym>EPSG:7030</ellipsoidacronym><geographicflag>false</"
    "geographicflag></spatialrefsys>";

static std::string opt(const char* k, const std::string& v) {
  return std::string("<Option name=\"") + k + "\" type=\"QString\" value=\"" +
         xml(v) + "\"/>";
}
static std::string prop(const char* k, const std::string& v) {
  return std::string("<prop k=\"") + k + "\" v=\"" + xml(v) + "\"/>";
}

static std::string marker_symbol(const std::string& name,
                                 const std::string& color, double size,
                                 const std::string& shape) {
  std::string kv[][2] = {
      {"color", color + ",255"},         {"name", shape},
      {"size", std::to_string(size)},    {"size_unit", "MM"},
      {"outline_color", "35,35,35,255"}, {"outline_width", "0.2"},
      {"outline_style", "solid"}};
  std::string s = "<symbol type=\"marker\" name=\"" + name +
                  "\" alpha=\"1\" clip_to_extent=\"1\" force_rhr=\"0\"><layer "
                  "class=\"SimpleMarker\" enabled=\"1\" "
                  "pass=\"0\" locked=\"0\"><Option type=\"Map\">";
  for (auto& p : kv)
    s += opt(p[0].c_str(), p[1]);
  s += "</Option>";
  for (auto& p : kv)
    s += prop(p[0].c_str(), p[1]);
  return s + "</layer></symbol>";
}

static std::string line_symbol(const std::string& name,
                               const std::string& color, double width) {
  std::string kv[][2] = {
      {"line_color", color + ",255"}, {"line_width", std::to_string(width)},
      {"line_width_unit", "MM"},      {"line_style", "solid"},
      {"capstyle", "round"},          {"joinstyle", "round"}};
  std::string s = "<symbol type=\"line\" name=\"" + name +
                  "\" alpha=\"0.9\" clip_to_extent=\"1\" "
                  "force_rhr=\"0\"><layer class=\"SimpleLine\" enabled=\"1\" "
                  "pass=\"0\" locked=\"0\"><Option type=\"Map\">";
  for (auto& p : kv)
    s += opt(p[0].c_str(), p[1]);
  s += "</Option>";
  for (auto& p : kv)
    s += prop(p[0].c_str(), p[1]);
  return s + "</layer></symbol>";
}

static std::string ramp_color(double t) {
  RGB c = colormap((float)(0.15 + 0.85 * t));
  return std::to_string(c.r) + "," + std::to_string(c.g) + "," +
         std::to_string(c.b);
}

static std::string renderer_xml(const LayerSpec& L) {
  std::ostringstream o;
  if (L.kind == LayerSpec::Lines) {
    o << "<renderer-v2 type=\"singleSymbol\" symbollevels=\"0\" "
         "enableorderby=\"0\" forceraster=\"0\"><symbols>"
      << line_symbol("0", L.color, L.size) << "</symbols></renderer-v2>";
  } else if (L.style == LayerSpec::Heatmap) {
    std::string stops = "0.25;" + ramp_color(0.25) + ",200:0.5;" +
                        ramp_color(0.5) + ",220:0.75;" + ramp_color(0.75) +
                        ",240";
    o << "<renderer-v2 type=\"heatmapRenderer\" radius=\"10\" "
         "radius_unit=\"0\" quality=\"3\" max_value=\"0\" "
         "weight_expression=\"&quot;"
      << xml(L.field)
      << "&quot;\" enableorderby=\"0\" forceraster=\"0\" "
         "radius_map_unit_scale=\"3x:0,0,0,0,0,0\">"
         "<colorramp type=\"gradient\" name=\"[source]\"><Option type=\"Map\">"
      << opt("color1", "0,0,4,0") << opt("color2", "252,255,164,255")
      << opt("discrete", "0") << opt("rampType", "gradient")
      << opt("stops", stops) << "</Option>" << prop("color1", "0,0,4,0")
      << prop("color2", "252,255,164,255") << prop("discrete", "0")
      << prop("rampType", "gradient") << prop("stops", stops)
      << "</colorramp></renderer-v2>";
  } else if (L.style == LayerSpec::Graduated && L.breaks.size() >= 2) {
    size_t n = L.breaks.size() - 1;
    o << "<renderer-v2 type=\"graduatedSymbol\" attr=\"" << xml(L.field)
      << "\" graduatedMethod=\"GraduatedColor\" symbollevels=\"0\" "
         "enableorderby=\"0\" forceraster=\"0\"><ranges>";
    for (size_t i = 0; i < n; ++i) {
      char b[256];
      std::snprintf(b, sizeof b,
                    "<range lower=\"%.6f\" upper=\"%.6f\" symbol=\"%zu\" "
                    "label=\"%.3f - %.3f\" render=\"true\"/>",
                    L.breaks[i], L.breaks[i + 1], i, L.breaks[i],
                    L.breaks[i + 1]);
      o << b;
    }
    o << "</ranges><symbols>";
    for (size_t i = 0; i < n; ++i)
      o << marker_symbol(std::to_string(i),
                         ramp_color(n > 1 ? (double)i / (n - 1) : 1.0), L.size,
                         L.marker);
    o << "</symbols><mode name=\"quantile\"/></renderer-v2>";
  } else {
    o << "<renderer-v2 type=\"singleSymbol\" symbollevels=\"0\" "
         "enableorderby=\"0\" forceraster=\"0\"><symbols>"
      << marker_symbol("0", L.color, L.size, L.marker)
      << "</symbols></renderer-v2>";
  }
  return o.str();
}

void write_qgs_project(const std::string& path, const ProjectSpec& p) {
  std::ofstream o(path);
  if (!o) die("cannot write %s", path.c_str());
  auto lid = [](size_t i) { return "planeloc_layer_" + std::to_string(i); };
  // EPSG:3857 extent for the canvas
  auto merc = [](double lon, double lat, double& x, double& y) {
    x = lon * kDeg * 6378137.0;
    y = std::log(std::tan(kPi / 4 + lat * kDeg / 2)) * 6378137.0;
  };
  double x0, y0, x1, y1;
  merc(p.minlon, p.minlat, x0, y0);
  merc(p.maxlon, p.maxlat, x1, y1);
  o << "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>\n"
    << "<qgis projectname=\"" << xml(p.title)
    << "\" version=\"3.34.0-Prizren\">\n"
    << "  <homePath path=\"\"/>\n  <title>" << xml(p.title) << "</title>\n"
    << "  <projectCrs>" << kSrs3857 << "</projectCrs>\n"
    << "  <layer-tree-group>\n    <customproperties/>\n";
  for (size_t i = 0; i < p.layers.size(); ++i) {
    const auto& L = p.layers[i];
    o << "    <layer-tree-layer id=\"" << lid(i) << "\" name=\"" << xml(L.name)
      << "\" source=\"./" << xml(L.file) << "\" providerKey=\"ogr\" checked=\""
      << (L.visible ? "Qt::Checked" : "Qt::Unchecked")
      << "\" expanded=\"0\"><customproperties/></layer-tree-layer>\n";
  }
  if (p.osm_basemap)
    o << "    <layer-tree-layer id=\"planeloc_osm\" name=\"OpenStreetMap\" "
         "source=\"\" providerKey=\"wms\" "
         "checked=\"Qt::Checked\" "
         "expanded=\"0\"><customproperties/></layer-tree-layer>\n";
  o << "  </layer-tree-group>\n  <projectlayers>\n";
  for (size_t i = 0; i < p.layers.size(); ++i) {
    const auto& L = p.layers[i];
    const char* g = L.kind == LayerSpec::Lines ? "Line" : "Point";
    const char* wkb = L.kind == LayerSpec::Lines ? "LineString" : "Point";
    o << "    <maplayer type=\"vector\" geometry=\"" << g << "\" wkbType=\""
      << wkb
      << "\" hasScaleBasedVisibilityFlag=\"0\" autoRefreshEnabled=\"0\" "
         "refreshOnNotifyEnabled=\"0\">\n"
      << "      <id>" << lid(i) << "</id>\n      <datasource>./" << xml(L.file)
      << "</datasource>\n"
      << "      <layername>" << xml(L.name) << "</layername>\n      <srs>"
      << kSrs4326 << "</srs>\n"
      << "      <provider encoding=\"UTF-8\">ogr</provider>\n      "
      << renderer_xml(L) << "\n"
      << "      "
         "<blendMode>0</blendMode><featureBlendMode>0</"
         "featureBlendMode><layerOpacity>1</layerOpacity>\n"
      << "    </maplayer>\n";
  }
  if (p.osm_basemap)
    o << "    <maplayer type=\"raster\" hasScaleBasedVisibilityFlag=\"0\" "
         "autoRefreshEnabled=\"0\">\n"
      << "      <id>planeloc_osm</id>\n"
      << "      "
         "<datasource>type=xyz&amp;url=https://tile.openstreetmap.org/%7Bz%7D/"
         "%7Bx%7D/%7By%7D.png&amp;zmax=19&amp;zmin=0</datasource>\n"
      << "      <layername>OpenStreetMap</layername>\n      <srs>" << kSrs3857
      << "</srs>\n"
      << "      <provider>wms</provider>\n    </maplayer>\n";
  o << "  </projectlayers>\n  <layerorder>\n";
  for (size_t i = 0; i < p.layers.size(); ++i)
    o << "    <layer id=\"" << lid(i) << "\"/>\n";
  if (p.osm_basemap) o << "    <layer id=\"planeloc_osm\"/>\n";
  char ext[256];
  std::snprintf(
      ext, sizeof ext,
      "<xmin>%.3f</xmin><ymin>%.3f</ymin><xmax>%.3f</xmax><ymax>%.3f</ymax>",
      x0, y0, x1, y1);
  o << "  </layerorder>\n  <mapcanvas name=\"theMapCanvas\" "
       "annotationsVisible=\"1\">\n    <units>meters</units>\n"
    << "    <extent>" << ext << "</extent>\n    <rotation>0</rotation>\n"
    << "    <destinationsrs>" << kSrs3857
    << "</destinationsrs>\n  </mapcanvas>\n"
    << "  <properties>\n    <Paths><Absolute "
       "type=\"bool\">false</Absolute></Paths>\n"
    << "    <Measurement><DistanceUnits "
       "type=\"QString\">meters</DistanceUnits></Measurement>\n  "
       "</properties>\n"
    << "</qgis>\n";
}

// --------------------------------------------------- spectrogram PNG ---
Image render_spectrogram(const Spectrogram& w, double t0, double t1, double f0,
                         double f1, const std::vector<Curve>& curves, int ys) {
  int fa = std::clamp((int)std::floor((t0 - w.t0) / w.hop), 0,
                      std::max(0, w.nframes - 1));
  int fb = std::clamp((int)std::ceil((t1 - w.t0) / w.hop), 0,
                      std::max(0, w.nframes - 1));
  int ba = std::clamp((int)std::floor(f0 / w.df) - w.bin0, 0,
                      std::max(0, w.nbins - 1));
  int bb = std::clamp((int)std::ceil(f1 / w.df) - w.bin0, 0,
                      std::max(0, w.nbins - 1));
  int W = std::max(1, fb - fa + 1), H = std::max(1, (bb - ba + 1) * ys);
  Image img(W, H);
  float mx = 1e-6f;
  for (int f = fa; f <= fb; ++f)
    for (int b = ba; b <= bb; ++b)
      mx = std::max(mx, w.at(f, b));
  mx = std::min(mx, 6.f);
  for (int x = 0; x < W; ++x)
    for (int y = 0; y < H; ++y) {
      int b = bb - y / ys;
      img.set(x, y, colormap(w.at(fa + x, b) / mx));
    }
  for (const auto& c : curves) {
    int px = -1, py = -1;
    for (auto& [t, f] : c.tf) {
      int x = (int)std::lround((t - w.t0) / w.hop) - fa;
      int y = (int)std::lround((bb - (f / w.df - w.bin0)) * ys);
      if (px >= 0 && std::abs(x - px) < 20) img.line(px, py, x, y, c.color);
      px = x;
      py = y;
    }
  }
  return img;
}

}  // namespace pl
