#pragma once

#include <string>
#include <vector>

#include <geometry.hh>
#include <harmonic.h>

namespace CGB {

using namespace Geometry;

class ConcaveGB {
  ConcaveGB();
  ~ConcaveGB();
public:
  void loadControlPoints(const std::string &filename);
  void loadControlPoints(size_t n, size_t d, std::vector<Point3D>);
  void generateDomain();
  TriMesh evaluate(double resolution) const;
private:
  Point3D evaluate(const Point2D &uv) const;

  using Ribbon = std::vector<std::vector<Point3D>>;
  size_t n_, d_;
  std::vector<Point2D> domain_;
  std::vector<HarmonicMap *> parameters_;
  std::vector<Ribbon> ribbons_;
};

}
