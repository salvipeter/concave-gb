#pragma once

#include <string>
#include <vector>

#include <geometry.hh>

struct HarmonicMap;

namespace CGB {

using namespace Geometry;

class ConcaveGB {
public:
  enum class CentralWeight { ORIGINAL, ZERO };
  using Ribbon = std::vector<PointVector>; // (degree + 1) * layer

  ConcaveGB();
  ~ConcaveGB();

  void setParamLevels(unsigned int levels);
  void setCentralWeight(CentralWeight type);
  void setCentralControlPoint(const Point3D &p);
  bool setRibbons(const std::vector<Ribbon> &ribbons, bool generate_domain = true);
  bool loadControlPoints(const std::string &filename, bool generate_domain = true);
  void generateDomain();
  TriMesh evaluate(double resolution) const;
private:
  DoubleVector localCoordinates(const Point2D &uv) const;
  double weight(const DoubleVector &bc, size_t i, size_t j, size_t k) const;
  Point3D evaluate(const DoubleVector &bc) const;

  unsigned int param_levels_;
  CentralWeight central_weight_;
  Point3D central_cp_;
  Point2DVector domain_;
  std::vector<HarmonicMap *> parameters_;
  std::vector<Ribbon> ribbons_;

  // Caching
  mutable double last_resolution_;
  mutable std::vector<DoubleVector> param_cache_;
  mutable TriMesh mesh_cache_;
};

}
