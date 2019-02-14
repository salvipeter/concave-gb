#pragma once

#include "cgb.hh"

namespace CGB {

class Harmonic;

class CurvedDomain {
public:
  CurvedDomain(const std::vector<ConcaveGB::Ribbon> &ribbons);
  std::unique_ptr<Harmonic> init(size_t levels) const;
  Point2DVector polyline(size_t resolution) const;
private:
  std::vector<Point2DVector> plane_curves_;
};

}
