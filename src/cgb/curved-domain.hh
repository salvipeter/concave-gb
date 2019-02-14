#pragma once

#include "cgb.hh"

namespace CGB {

class Harmonic;

class CurvedDomain {
public:
  CurvedDomain(const std::vector<ConcaveGB::Ribbon> &ribbons, size_t levels);
  Point2DVector polyline(size_t resolution) const;
  DoubleVector harmonic(const Point2D &p) const;
private:
  std::vector<Point2DVector> plane_curves_;
  std::unique_ptr<Harmonic> harmonic_;
};

}
