// Optional native QGIS integration: links libqgis_core and writes a .qgz project through the
// QGIS C++ API (QgsProject / QgsVectorLayer / renderers). Enabled with -DPLANELOC_WITH_QGIS=ON.
#include "export.h"

#include <cmath>
#include <cstdio>

#ifdef PLANELOC_WITH_QGIS
#include <qgsapplication.h>
#include <qgscoordinatereferencesystem.h>
#include <qgscolorramp.h>
#include <qgsgraduatedsymbolrenderer.h>
#include <qgsheatmaprenderer.h>
#include <qgslayertree.h>
#include <qgslayertreelayer.h>
#include <qgslinesymbol.h>
#include <qgsmarkersymbol.h>
#include <qgsproject.h>
#include <qgsprojectviewsettings.h>
#include <qgsreferencedgeometry.h>
#include <qgsrasterlayer.h>
#include <qgsrendererrange.h>
#include <qgssinglesymbolrenderer.h>
#include <qgsvectorlayer.h>

#include <QColor>
#include <QString>
#include <QVariantMap>
#endif

namespace pl {

#ifdef PLANELOC_WITH_QGIS

bool qgis_native_available() { return true; }

static QColor qcolor(const std::string& rgb) {
  int r = 255, g = 0, b = 0;
  std::sscanf(rgb.c_str(), "%d,%d,%d", &r, &g, &b);
  return QColor(r, g, b);
}

bool write_qgis_native(const std::string& path, const std::string& dir, const ProjectSpec& p) {
  static int argc = 1;
  static char arg0[] = "planeloc";
  static char* argv[] = {arg0, nullptr};
  QgsApplication app(argc, argv, false);
  QgsApplication::initQgis();
  {
    QgsProject project;
    project.setTitle(QString::fromStdString(p.title));
    project.setCrs(QgsCoordinateReferenceSystem(QStringLiteral("EPSG:3857")));
    QList<QgsMapLayer*> layers;
    for (const LayerSpec& L : p.layers) {
      QString src = QString::fromStdString(dir + "/" + L.file);
      auto* vl = new QgsVectorLayer(src, QString::fromStdString(L.name), QStringLiteral("ogr"));
      if (!vl->isValid()) {
        delete vl;
        continue;
      }
      if (L.kind == LayerSpec::Lines) {
        QVariantMap props;
        props[QStringLiteral("line_color")] = QString::fromStdString(L.color);
        props[QStringLiteral("line_width")] = QString::number(L.size);
        vl->setRenderer(new QgsSingleSymbolRenderer(QgsLineSymbol::createSimple(props)));
      } else if (L.style == LayerSpec::Heatmap) {
        auto* hm = new QgsHeatmapRenderer();
        hm->setWeightExpression(QStringLiteral("\"%1\"").arg(QString::fromStdString(L.field)));
        hm->setRadius(10);
        hm->setColorRamp(new QgsGradientColorRamp(QColor(0, 0, 4, 0), QColor(252, 255, 164)));
        vl->setRenderer(hm);
      } else if (L.style == LayerSpec::Graduated && L.breaks.size() >= 2) {
        QgsRangeList ranges;
        size_t n = L.breaks.size() - 1;
        for (size_t i = 0; i < n; ++i) {
          RGB c = colormap((float)(0.15 + 0.85 * (n > 1 ? (double)i / (n - 1) : 1.0)));
          QVariantMap props;
          props[QStringLiteral("color")] = QStringLiteral("%1,%2,%3").arg(c.r).arg(c.g).arg(c.b);
          props[QStringLiteral("size")] = QString::number(L.size);
          props[QStringLiteral("name")] = QString::fromStdString(L.marker);
          ranges.append(QgsRendererRange(L.breaks[i], L.breaks[i + 1], QgsMarkerSymbol::createSimple(props),
                                         QStringLiteral("%1 - %2").arg(L.breaks[i], 0, 'f', 3).arg(L.breaks[i + 1], 0, 'f', 3)));
        }
        vl->setRenderer(new QgsGraduatedSymbolRenderer(QString::fromStdString(L.field), ranges));
      } else {
        QVariantMap props;
        props[QStringLiteral("color")] = QString::fromStdString(L.color);
        props[QStringLiteral("size")] = QString::number(L.size);
        props[QStringLiteral("name")] = QString::fromStdString(L.marker);
        vl->setRenderer(new QgsSingleSymbolRenderer(QgsMarkerSymbol::createSimple(props)));
      }
      layers.append(vl);
    }
    if (p.osm_basemap) {
      auto* osm = new QgsRasterLayer(
          QStringLiteral("type=xyz&url=https://tile.openstreetmap.org/%7Bz%7D/%7Bx%7D/%7By%7D.png&zmax=19&zmin=0"),
          QStringLiteral("OpenStreetMap"), QStringLiteral("wms"));
      if (osm->isValid()) layers.append(osm);
      else delete osm;
    }
    // Layers are inserted at the top of the layer tree, so add bottom-most (basemap) first.
    for (int i = layers.size() - 1; i >= 0; --i) project.addMapLayer(layers[i]);
    for (int i = 0; i < layers.size(); ++i) {
      bool vis = true;
      for (const LayerSpec& L : p.layers)
        if (QString::fromStdString(L.name) == layers[i]->name()) vis = L.visible;
      if (auto* node = project.layerTreeRoot()->findLayer(layers[i]->id())) node->setItemVisibilityChecked(vis);
    }
    auto merc = [](double lon, double lat, double& x, double& y) {
      x = lon * 0.017453292519943295 * 6378137.0;
      y = std::log(std::tan(0.7853981633974483 + lat * 0.008726646259971648)) * 6378137.0;
    };
    double x0, y0, x1, y1;
    merc(p.minlon, p.minlat, x0, y0);
    merc(p.maxlon, p.maxlat, x1, y1);
    project.viewSettings()->setDefaultViewExtent(
        QgsReferencedRectangle(QgsRectangle(x0, y0, x1, y1), QgsCoordinateReferenceSystem(QStringLiteral("EPSG:3857"))));
    bool ok = project.write(QString::fromStdString(path));
    project.clear();
    QgsApplication::exitQgis();
    return ok;
  }
}

#else

bool qgis_native_available() { return false; }
bool write_qgis_native(const std::string&, const std::string&, const ProjectSpec&) { return false; }

#endif

}  // namespace pl
