#pragma once

#include <string>
#include <vector>

#include <geometry.hh>

struct HarmonicMap;

namespace CGB {

using namespace Geometry;

class ConcaveGB {
public:
  enum class CentralWeight { ORIGINAL = 0, ZERO = 1, NTH = 2, HARMONIC = 3 };
  using Ribbon = std::vector<PointVector>; // (degree + 1) * layer

  ConcaveGB();
  ~ConcaveGB();

  // Level of detail for the discrete approximation of harmonic coordinates.
  // The program uses a set of n bitmaps, each of dimension 2^levels x 2^levels.
  // Usual values are between 9 and 11.
  void setParamLevels(unsigned int levels);

  // Sets the type of central weight computation.
  // ORIGINAL: the total weight deficiency (as described in the original paper)
  // ZERO: zero, using normalization of the weights
  // NTH: 1/n-th of the total weight deficiency (and all weights are normalized)
  // HARMONIC: computed from harmonic coordinates (and all weights are normalized)
  void setCentralWeight(CentralWeight type);

  // The given tolerance controls the minimum distance between non-adjacent segments.
  // Note that the domain is always normalized to the [-1,1]x[-1,1] square.
  void setDomainTolerance(double tolerance);

  // Alternative parameterization - positive values result in smaller weight deficiency.
  void setParameterDilation(double dilation);

  // Sets the central control point position.
  void setCentralControlPoint(const Point3D &p);

  // Sets the ribbons. Note that this does not set the central control point.
  // If generate_domain is set to false, the program assumes that the domain can stay the same.
  // A Ribbon consists of rows of control points, each row contains degree + 1 control points.
  bool setRibbons(const std::vector<Ribbon> &ribbons, bool generate_domain = true);

  // Reads model-specific options. The file format is the following:
  //   central_weight (as an integer)
  //   domain_tolerance
  //   parameter_dilation
  bool loadOptions(std::istream &is);

  // Sets the ribbons and the central control point.
  // If generate_domain is set to false, the program assumes that the domain can stay the same.
  // The file format is the following:
  //   n
  //   c[x] c[y] c[z]
  //   d_1 l_1
  //   p_101[x] p_101[y] p101[z]
  //   ...
  //   d_2 l_2
  //   p_201[x] p_201[y] p_201[z]
  //   ...
  // where n is the # of ribbons, c is the central control point,
  //       d_i and l_i are the degree and layer of the i-th ribbon, respectively,
  //       p_ijk is the control point of the i-th ribbon in the j-th column of the k-th row
  // Points are stored row by row, so
  //   (1,0,1)->...->(1,d1,1)->(1,0,2)->...->(1,d1,l1)->(2,0,1)->...->(n,dn,ln)
  // Here i and k is indexed from 1, while j is indexed from 0, for convenience.
  bool loadControlPoints(std::istream &is, bool generate_domain = true);

  // Regenerates the domain and updates the harmonic maps.
  // Assumes that the ribbons are set.
  void generateDomain();

  // Evaluates the surface at a given resolution, resulting in a triangle mesh.
  // If the resolution did not change since the last call,
  // cached parameters and mesh topology are used.
  // If the resolution is positive, it gives the maximal area of a triangle in the domain.
  // Since the domain is always in [-1,1]x[-1,1], usual values range betwen 1e-5 and 1e-3.
  // If the resolution is negative, a bitmap-based triangulation is used,
  // where the resolution means downsampling from the ParamLevels resolution.
  // So e.g. ParamLevels = 9 (the default) is a 2^9=512x512 bitmap,
  // with resolution = -2 the triangulation will use a 2^(9-2)=2^7=128x128 bitmap.
  TriMesh evaluate(double resolution) const;

private:
  Point2D barycentricSD(const DoubleVector &bc, size_t i) const;
  DoubleVector localCoordinates(const Point2D &uv) const;
  double weight(const DoubleVector &bc, size_t i, size_t j, size_t k) const;
  Point3D evaluate(const DoubleVector &bc) const;

  unsigned int param_levels_;
  CentralWeight central_weight_;
  double domain_tolerance_;
  double parameter_dilation_;
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
