#pragma once

#include <string>
#include <vector>

#include <geometry.hh>

struct mec_t;

namespace CGB {

class CurvedDomain;
class Harmonic;

using namespace Geometry;

class ConcaveGB {
public:
  enum class CentralWeight { ORIGINAL = 0, ZERO = 1, NTH = 2, HARMONIC = 3 };
  enum class DomainType { PROJECTION = 0, NORMAL = 1, CORRECTED = 2, CURVED = 3 };
  using Ribbon = std::vector<PointVector>; // (degree + 1) * layer

  ConcaveGB();
  ~ConcaveGB();

  // Level of detail for the discrete approximation of harmonic coordinates.
  // The program uses a set of n bitmaps, each of dimension 2^levels x 2^levels.
  // Usual values are between 9 and 11.
  void setParamLevels(size_t levels);

  // Sets the type of central weight computation.
  // ORIGINAL: the total weight deficiency (as described in the original paper)
  // ZERO: zero, using normalization of the weights
  // NTH: 1/n-th of the total weight deficiency (and all weights are normalized)
  // HARMONIC: computed from harmonic coordinates (and all weights are normalized)
  void setCentralWeight(CentralWeight type);

  // Sets the weight of the parameterization of concave corners (default = 1)
  void setConcaveWeight(double weight);

  // Sets the type of domain.
  // PROJECTION: a simple projection of the vertices into a best-fit plane
  // NORMAL: length-angle based (as in the paper)
  // CORRECTED: same, but using outer angles for the obtuse angles (when possible)
  // CURVED: curved domain
  void setDomainType(DomainType type);

  // The given tolerance controls the minimum distance between non-adjacent segments.
  // Note that the domain is always normalized to the [0,1]x[0,1] square.
  void setDomainTolerance(double tolerance);

  // Alternative parameterization - positive values result in smaller weight deficiency.
  void setParameterDilation(double dilation);

  // Use maximum entropy coordinates (false: use harmonic coordinates).
  void setMaxEnt(bool maxent);

  // Sets the central control point position.
  void setCentralControlPoint(const Point3D &p);

  // Add a dummy control point at all concave corners
  void setFillConcaveCorners(bool fill_concave);

  // Sets the ribbons. Note that this does not set the central control point.
  // If generate_domain is set to false, the program assumes that the domain can stay the same.
  // A Ribbon consists of rows of control points, each row contains degree + 1 control points.
  bool setRibbons(const std::vector<Ribbon> &ribbons, bool generate_domain = true);

  // Reads model-specific options. The file format is the following:
  //   central_weight (as an integer)
  //   domain_type
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
  // Resolution gives the maximal area of a triangle in the domain.
  // Since the domain is always in [0,1]x[0,1], usual values range betwen 1e-5 and 1e-3.
  TriMesh evaluate(double resolution) const;

  // Same as evaluate(), but a bitmap-based triangulation is used,
  // where the triangulation uses a 2^resolution x 2^resolution bitmap.
  TriMesh evaluateRegular(size_t resolution) const;

private:
  Point2DVector generateSimilarityDomain() const;
  Point2DVector generateProjectedDomain() const;
  void generateDelaunayMesh(double resolution) const;
  void generateRegularMesh(size_t downsampling) const;
  DoubleVector localCoordinates(const Point2D &uv) const;
  Point3D evaluate(const DoubleVector &bc) const;
  TriMesh evaluateImpl() const;

  size_t param_levels_;
  CentralWeight central_weight_;
  double concave_weight_;
  DomainType domain_type_;
  double domain_tolerance_;
  double parameter_dilation_;
  bool use_maxent_;
  bool fill_concave_corners_;
  Point3D central_cp_;
  Point2DVector domain_;
  std::unique_ptr<CurvedDomain> curved_domain_;
  std::vector<bool> concave_;
  std::unique_ptr<Harmonic> parameters_;
  mec_t *mec_parameters_;
  std::vector<Ribbon> ribbons_;

  struct ExtraCP {
    size_t i;
    Point3D p;
  };
  std::vector<ExtraCP> extra_cp_;

  // Just for testing
  mutable double def;
  mutable PointVector def_points;

  // Caching
  mutable double last_mesh_size_;
  mutable bool last_regular_;
  mutable std::vector<DoubleVector> param_cache_;
  mutable TriMesh mesh_cache_;
};

}
