#include "bezier.hh"
#include "curved-domain.hh"
#include "harmonic.hh"
#include "lsq-plane.hh"

namespace CGB {

static void scalePoints(Point2DVector &points) {
  Point2D min = points[0], max = points[0];
  for (size_t i = 1; i < points.size(); ++i) {
    for (size_t j = 0; j < 2; ++j) {
      min[j] = std::min(min[j], points[i][j]);
      max[j] = std::max(max[j], points[i][j]);
    }
  }
  auto d = max - min;
  double margin = 0.025;
  double len = std::max(d[0], d[1]) * (1.0 + 2.0 * margin);
  for (auto &p : points) {
    p[0] = (p[0] - min[0]) / len + margin;
    p[1] = (p[1] - min[1]) / len + margin;
  }
}

CurvedDomain::CurvedDomain(const std::vector<ConcaveGB::Ribbon> &ribbons, size_t levels) {
  PointVector pv;
  for (const auto &r : ribbons) {
    const auto &cp = r.front();
    for (size_t i = 1; i < cp.size(); ++i)
      pv.push_back(cp[i]);
  }
  auto projected = LSQPlane::projectToBestFitPlane(pv);
  scalePoints(projected);

  int index = -1;
  for (const auto &r : ribbons) {
    size_t d = r.front().size() - 1;
    Point2DVector cp;
    cp.push_back(index >= 0 ? projected[index] : projected.back());
    for (size_t i = 1; i <= d; ++i)
      cp.push_back(projected[++index]);
    plane_curves_.push_back(cp);
  }

  harmonic_ = std::make_unique<Harmonic>(plane_curves_, levels);
}

Point2DVector
CurvedDomain::polyline(size_t resolution) const {
  Point2DVector result;
  for (const auto &c : plane_curves_)
    for (size_t i = 0; i < resolution; ++i) {
      double u = (double)i / resolution;
      result.push_back(bezierEval(c, u));
    }
  return result;
}

DoubleVector
CurvedDomain::harmonic(const Point2D &p) const {
  return harmonic_->eval(p);
}

}
