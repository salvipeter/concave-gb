#pragma once

#include <geometry.hh>

namespace CGB {

using namespace Geometry;

struct GridValue {
  bool boundary;
  double value;
};
using HarmonicMap = std::vector<GridValue>;

class Harmonic {
public:
  Harmonic(const std::vector<Point2DVector> &curves, size_t levels);
  DoubleVector eval(const Point2D &uv) const;
private:
  size_t levels_, size_;
  std::vector<HarmonicMap> maps_;
};

}
