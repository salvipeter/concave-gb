#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stack>

#include <harmonic.h>

#include "lsq-plane.hh"

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include <triangle.h>
}

#include "cgb.hh"

#ifndef M_PI
static const double M_PI = std::acos(-1);
#endif

static const double EPSILON = 1.0e-5;

namespace CGB {

ConcaveGB::ConcaveGB() :
  param_levels_(9), central_weight_(CentralWeight::ZERO), domain_tolerance_(0.2),
  parameter_dilation_(0.0), central_cp_(0, 0, 0), last_resolution_(std::nan(""))
{
}

ConcaveGB::~ConcaveGB() {
  for (auto p : parameters_)
    harmonic_free(p);
}


void
ConcaveGB::setParamLevels(size_t levels) {
  param_levels_ = levels;
}

void
ConcaveGB::setCentralWeight(ConcaveGB::CentralWeight type) {
  central_weight_ = type;
}

void
ConcaveGB::setDomainTolerance(double tolerance) {
  domain_tolerance_ = tolerance;
}

void
ConcaveGB::setParameterDilation(double dilation) {
  parameter_dilation_ = dilation;
}

void
ConcaveGB::setCentralControlPoint(const Point3D &p) {
  central_cp_ = p;
}

bool
ConcaveGB::setRibbons(const std::vector<Ribbon> &ribbons, bool generate_domain) {
  ribbons_ = ribbons;

  if (generate_domain)
    generateDomain();

  return true;
}

bool
ConcaveGB::loadOptions(std::istream &is) {
  unsigned int cw;
  is >> cw >> domain_tolerance_ >> parameter_dilation_;

  if (cw > 3 || domain_tolerance_ < 0 || parameter_dilation_ < 0)
    return false;

  central_weight_ = static_cast<CentralWeight>(cw);

  return is.good();
}

bool
ConcaveGB::loadControlPoints(std::istream &is, bool generate_domain) {
  ribbons_.clear();

  size_t n, d, l;

  is >> n;
  Point3D p;
  is >> central_cp_[0] >> central_cp_[1] >> central_cp_[2];
  ribbons_.reserve(n);
  for (size_t side = 0; side < n; ++side) {
    is >> d >> l;
    Ribbon r; r.reserve(l);
    for (size_t row = 0; row < l; ++row) {
      PointVector one_row; one_row.reserve(d + 1);
      for (size_t col = 0; col <= d; ++col) {
        is >> p[0] >> p[1] >> p[2];
        one_row.push_back(p);
      }
      r.push_back(std::move(one_row));
    }
    ribbons_.push_back(std::move(r));
  }
  
  if (generate_domain && is.good())
    generateDomain();

  return is.good();
}

namespace {

  // Forces x into [min, max].
  double
  inrange(double min, double x, double max) {
    if (x < min)
      return min;
    if (x > max)
      return max;
    return x;
  }

  // Normalize angles, but consider concave angles as the negative of the outer (smaller) angles.
  // This can fail if the sum is close to 0, or if the target is (exactly) 0.
  bool
  normalizeSmallerAngles(DoubleVector &angles) {
    auto flip = [](double alpha) { return alpha < 0 ? -(M_PI + alpha) : M_PI - alpha; };
    size_t concave_count = 0;
    double angle_sum = 0.0;
    for (const auto &angle : angles) {
      if (angle < 0)
        ++concave_count;
      angle_sum += flip(angle);
    }
    double target = (angles.size() - 2 - 2 * concave_count) * M_PI;
    if (target == 0 || std::abs(angle_sum) < EPSILON)
      return false;
    double angle_multiplier = target / angle_sum;
    for (auto &angle : angles)
      angle = flip(flip(angle) * angle_multiplier);
    return true;
  }

  // Normalize angles the old-fashioned way.
  // This also enlarges concave angles, which is somewhat counter-intuitive,
  // but we need it as a fallback when normalizeSmallerAngles fails.
  void
  normalizeInnerAngles(DoubleVector &angles) {
    auto flip = [](double alpha) { return M_PI - alpha; };
    double angle_sum = 0.0;
    for (const auto &angle : angles)
      angle_sum += flip(angle);
    double target = (angles.size() - 2) * M_PI;
    double angle_multiplier = target / angle_sum;
    for (auto &angle : angles)
      angle = flip(flip(angle) * angle_multiplier);
  }

  // Rescale to [-1,-1]x[1,1].
  void
  rescaleDomain(Point2DVector &domain) {
    double minx = 0.0, miny = 0.0, maxx = 0.0, maxy = 0.0;
    for (const auto &v : domain) {
      minx = std::min(minx, v[0]); miny = std::min(miny, v[1]);
      maxx = std::max(maxx, v[0]); maxy = std::max(maxy, v[1]);
    }
    double width = std::max(maxx - minx, maxy - miny);
    Point2D topleft(-1.0, -1.0), minp(minx, miny);
    for (auto &v : domain)
      v = topleft + (v - minp) * 2.0 / width;
  }

  // Generates a domain in [-1,1]x[-1,1], but does not check for validity.
  Point2DVector
  generateAngleLengthDomain(const DoubleVector &angles, const DoubleVector &lengths) {
    size_t n = angles.size();

    // Compute open domain
    Point2DVector domain(n);
    double dir = 0.0;
    Vector2D prev_v(0.0, 0.0);
    for (size_t i = 0; i < n; ++i) {
      domain[i] = prev_v + Vector2D(std::cos(dir), std::sin(dir)) * lengths[i];
      dir += angles[i];
      prev_v = domain[i];
    }

    // Compute closed domain
    double length_sum = std::accumulate(lengths.begin(), lengths.end(), 0.0);
    double accumulated = 0.0;
    for (size_t i = 0; i < n; ++i) {
      accumulated += lengths[i];
      domain[i] -= domain.back() * accumulated / length_sum;
    }

    rescaleDomain(domain);

    return domain;
  }

  using Segment = Point2D[2];

  double segmentSegmentDistance(const Segment &s1, const Segment &s2) {
    const auto &p1 = s1[0], &p2 = s2[0];
    auto v1 = s1[1] - p1, v2 = s2[1] - p2;
    auto d = p1 - p2;
    auto normal = [](const Vector2D &v) { return Vector2D(-v[1], v[0]).normalize(); };
    double denom = v2[0] * v1[1] - v2[1] * v1[0];

    if (std::abs(denom) < EPSILON) {
      // Parallel
      if (v1 * v2 < 0 && v1 * d >= 0)
        return d.norm();
      if (v1 * v2 > 0 && v1 * (d + v1) <= 0)
        return (d + v1).norm();
      if (v1 * v2 > 0 && v1 * (d - v2) >= 0)
        return (d - v2).norm();
      if (v1 * v2 < 0 && v1 * (d + v1 - v2) <= 0)
        return (d + v1 - v2).norm();
      return std::abs(normal(v1) * d);
    }

    // Not parallel
    double t1 = (v2[1] * d[0] - v2[0] * d[1]) / denom;
    double t2 = std::abs(v2[0]) > std::abs(v2[1])
      ? (t1 * v1[0] + d[0]) / v2[0]
      : (t1 * v1[1] + d[1]) / v2[1];

    // Check for intersection
    if (0 < t1 && t1 < 1 && 0 < t2 && t2 < 1)
      return 0;

    // Try to find an endpoint close to an interval
    double u1 = -(d * v1) / v1.normSqr();          // s2[0]'s parameter when projected on s1
    if (0 < u1 && u1 < 1 && t2 <= 0)
      return std::abs(normal(v1) * d);
    double w1 = -((d - v2) * v1) / v1.normSqr();   // s2[1]'s parameter when projected on s1
    if (0 < w1 && w1 < 1 && t2 >= 1)
      return std::abs(normal(v1) * (d - v2));
    double u2 = (d * v2) / v2.normSqr();           // s1[0]'s parameter when projected on s2
    if (0 < u2 && u2 < 1 && t1 <= 0)
      return std::abs(normal(v2) * d);
    double w2 = ((d + v1) * v2) / v2.normSqr();    // s1[1]'s parameter when projected on s2
    if (0 < w2 && w2 < 1 && t1 >= 1)
      return std::abs(normal(v2) * (d + v1));

    // Find the best endpoint pair
    if (u1 <= 0 && u2 <= 0)
      return d.norm();
    if (w1 <= 0 && u2 >= 1)
      return (d - v2).norm();
    if (u1 >= 1 && w2 <= 0)
      return (d + v1).norm();
    if (w1 >= 1 && w2 >= 1)
      return (d + v1 - v2).norm();

    // Should not come here
    return 0;
  }

  // A domain is valid if the smallest distance between non-adjacent segments is above tolerance.
  bool
  isDomainValid(const Point2DVector &domain, double tolerance) {
    size_t n = domain.size();
    for (size_t i = 0; i < n - 2; ++i) {
      size_t im = (i + n - 1) % n;
      Segment si = { domain[im], domain[i] };
      for (size_t j = i + 2; j < n; ++j) {
        size_t jm = (j + n - 1) % n;
        Segment sj = { domain[jm], domain[j] };
        if ((i > 0 || j < n - 1) &&
            segmentSegmentDistance(si, sj) < tolerance)
          return false;
      }
    }
    return true;
  }

  // Multiplies all convex angles by a factor (1.1 for now),
  // and distributes the extra amount between the concave vertices uniformly.
  // If there are no concave angles, it returns false.
  bool
  enlargeDomainAngles(DoubleVector &angles) {
    const double factor = 1.1;
    auto flip = [](double alpha) { return M_PI - alpha; };
    size_t concave_count = 0;
    double angle_sum = 0.0;
    for (const auto &angle : angles) {
      if (angle < 0) {
        ++concave_count;
        angle_sum += flip(angle);
      } else
        angle_sum += flip(angle) * factor;
    }
    if (concave_count == 0)
      return false;
    double target = (angles.size() - 2) * M_PI;
    double delta = (target - angle_sum) / concave_count;
    for (auto &angle : angles)
      angle = angle < 0 ? flip(flip(angle) + delta) : flip(flip(angle) * factor);
    return true;
  }

}

// Generates a domain that mimics the angles and arc lengths of the boundaries.
Point2DVector
ConcaveGB::generateSimilarityDomain() const {
  Point2DVector domain;
  size_t n = ribbons_.size();

  // Compute lengths
  DoubleVector lengths; lengths.reserve(n);
  std::transform(ribbons_.begin(), ribbons_.end(), std::back_inserter(lengths),
                 [](const auto &r) { return BCurve(r[0]).arcLength(0.0, 1.0); });

  // Compute angles
  VectorVector der;
  DoubleVector angles; angles.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    size_t ip = (i + 1) % n;
    auto &r1 = ribbons_[i][0];
    auto &r2 = ribbons_[ip][0];
    size_t d = r1.size() - 1;
    auto v1 = (r1[d] - r1[d-1]).normalize();
    auto v2 = (r2[1] - r2[0]).normalize();
    double alpha = std::acos(inrange(-1, v1 * v2, 1));
    if (v1 * (ribbons_[ip][1][0] - r2[0]) > 0)
      alpha *= -1;
    angles.push_back(alpha);
  }

  if (!normalizeSmallerAngles(angles))
    normalizeInnerAngles(angles);

  domain = generateAngleLengthDomain(angles, lengths);

  while (!isDomainValid(domain, domain_tolerance_)) {
    enlargeDomainAngles(angles);
    domain = generateAngleLengthDomain(angles, lengths);
  }

  return domain;
}

// Generates a domain by projecting the corner vertices in a LSQ-fit plane.
Point2DVector
ConcaveGB::generateProjectedDomain() const {
  size_t n = ribbons_.size();
  PointVector pv; pv.reserve(n);
  for (const auto &r : ribbons_)
    pv.push_back(r[0].back());

  Point2DVector domain = LSQPlane::projectToBestFitPlane(pv);
  rescaleDomain(domain);

  return domain;
}

void
ConcaveGB::generateDomain() {
  // Invalidate everything
  domain_.clear();
  for (auto p : parameters_)
    harmonic_free(p);
  parameters_.clear();
  last_resolution_ = std::nan("");
  param_cache_.clear();
  mesh_cache_.clear();

  if (domain_tolerance_ > 0)
    domain_ = generateSimilarityDomain();
  else
    domain_ = generateProjectedDomain();

  {
    auto scale = [](Point2D p) {
      return (p + Point2D(1,1)) * 250 + Point2D(50,50);
    };
    std::ofstream f("domain.eps");
    f << "newpath\n";
    auto q = scale(domain_.back());
    f << q[0] << ' ' << q[1] << " moveto\n";
    for (const auto &p : domain_) {
      auto q = scale(p);
      f << q[0] << ' ' << q[1] << " lineto\n";
    }
    f << "stroke\nshowpage" << std::endl;
  }

  // Setup parameterization
  size_t n = ribbons_.size();
  DoubleVector points; points.reserve(3 * n);
  for (const auto &p : domain_) {
    points.push_back(p[0]);
    points.push_back(p[1]);
    points.push_back(0);
  }
  for (size_t i = 0; i < n; ++i) {
    points[3*i+2] = 1;
    auto map = harmonic_init(n, &points[0], param_levels_, EPSILON, false);
    parameters_.push_back(map);
    points[3*i+2] = 0;
  }
}

namespace {

  // B^n_i(u) for all i = 0..n
  void
  bernstein(size_t n, double u, DoubleVector &coeff) {
    coeff.clear(); coeff.reserve(n + 1);
    coeff.push_back(1.0);
    double u1 = 1.0 - u;
    for (size_t j = 1; j <= n; ++j) {
      double saved = 0.0;
      for (size_t k = 0; k < j; ++k) {
        double  tmp = coeff[k];
        coeff[k] = saved + tmp * u1;
        saved = tmp * u;
      }
      coeff.push_back(saved);
    }
  }

  // Returns the index of the segment closest to p.
  // An edge is defined by vertices i-1 and i.
  // (This computes the distances from _lines_, not segments.)
  size_t
  closestEdge(const Point2DVector &points, const Point2D &p) {
    size_t n = points.size();
    size_t result = 0;
    double min = -1;
    for (size_t i = 0; i < n; ++i) {
      size_t im = (i + n - 1) % n;
      Vector2D dev = p - points[im];
      Vector2D dir = (points[i] - points[im]).normalize();
      double dist = (dev - dir * (dev * dir)).norm();
      if (min < 0 || dist < min) {
        min = dist;
        result = i;
      }
    }
    return result;
  }

  // Returns the (s,d) system for side i, given the barycentric coordinates bc.
  Point2D
  barycentricSD(const DoubleVector &bc, size_t i, double dilation) {
    size_t n = bc.size(), im = (i + n - 1) % n;
    size_t imm = (i + n - 2) % n, ip = (i + 1) % n;
    double him = bc[im], hi = bc[i];
    double himm = bc[imm], hip = bc[ip];
    double d = (1.0 - him - hi) * (1.0 - dilation * himm * hip);
    double s = him + hi;
    if (std::abs(s) > EPSILON)
      s = hi / s;
    return Point2D(s, d);
  }

  // Returns mu^i_{j,k}, given the local coordinates sds (side i, row k, column j).
  double
  weight(size_t d, const Point2DVector &sds, size_t i, size_t j, size_t k) {
    size_t n = sds.size();
    double mu = 1.0;
    if (k < 2 && j < 2) {
      // alpha
      size_t im = (i + n - 1) % n;
      double di2 = std::pow(sds[i][1], 2), dim2 = std::pow(sds[im][1], 2);
      double denom = dim2 + di2;
      if (denom < EPSILON)
        mu = 0.5;
      else
        mu = dim2 / denom;
    } else if (k < 2 && j > d - 2) {
      // beta
      size_t ip = (i + 1) % n;
      double di2 = std::pow(sds[i][1], 2), dip2 = std::pow(sds[ip][1], 2);
      double denom = dip2 + di2;
      if (denom < EPSILON)
        mu = 0.5;
      else
        mu = dip2 / denom;
    } else if (j < k || j > d - k)
      mu = 0.0;
    else if (j == k || j == d - k)
      mu = 0.5;
    return mu;
  }

}

// Returns the barycentric coordinates for an (u,v) point in the domain.
DoubleVector
ConcaveGB::localCoordinates(const Point2D &uv) const {
  size_t n = parameters_.size();
  DoubleVector result(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    double value;
    if (!harmonic_eval(parameters_[i], const_cast<double *>(uv.data()), &value)) {
      // linear interpolation
      size_t k = closestEdge(domain_, uv), km = (k + n - 1) % n;
      double x = (uv - domain_[km]).norm() / (domain_[k] - domain_[km]).norm();
      result[km] = 1 - x;
      result[k] = x;
      return result;
    }
    result[i] = value;
  }
  return result;
}

// Generates a new mesh cache using Shewchuk's Triangle library.
void
ConcaveGB::generateDelaunayMesh(double resolution) const {
  param_cache_.clear();
  mesh_cache_.clear();

  // Input points
  size_t n = domain_.size();
  DoubleVector points; points.reserve(2 * n);
  for (const auto &p : domain_) {
    points.push_back(p[0]);
    points.push_back(p[1]);
  }

  // Input segments : just a closed polygon
  std::vector<int> segments; segments.reserve(2 * n);
  for (size_t i = 0; i < n; ++i) {
    segments[2*i]   = i;
    segments[2*i+1] = i + 1;
  }
  segments[2*n-1] = 0;
    
  // Setup output data structure
  struct triangulateio in, out;
  in.pointlist = &points[0];
  in.numberofpoints = n;
  in.numberofpointattributes = 0;
  in.pointmarkerlist = nullptr;
  in.segmentlist = &segments[0];
  in.numberofsegments = n;
  in.segmentmarkerlist = nullptr;
  in.numberofholes = 0;
  in.numberofregions = 0;

  // Setup output data structure
  out.pointlist = nullptr;
  out.pointattributelist = nullptr;
  out.pointmarkerlist = nullptr;
  out.trianglelist = nullptr;
  out.triangleattributelist = nullptr;
  out.segmentlist = nullptr;
  out.segmentmarkerlist = nullptr;

  // Call the library function [with maximum triangle area = resolution]
  std::ostringstream cmd;
  cmd << "pqa" << std::fixed << resolution << "DBPzQ";
  triangulate(const_cast<char *>(cmd.str().c_str()), &in, &out, (struct triangulateio *)nullptr);

  // Process the result
  param_cache_.reserve(out.numberofpoints);
  for (int i = 0; i < out.numberofpoints; ++i)
    param_cache_.emplace_back(localCoordinates(Point2D(out.pointlist[2*i],
                                                       out.pointlist[2*i+1])));
  mesh_cache_.resizePoints(param_cache_.size());
  for (int i = 0; i < out.numberoftriangles; ++i)
    mesh_cache_.addTriangle(out.trianglelist[3*i+0],
                            out.trianglelist[3*i+1],
                            out.trianglelist[3*i+2]);
}

namespace {

  void floodFill(std::vector<bool> &grid, size_t n, size_t x, size_t y) {
    using Coord = std::pair<size_t, size_t>;
    std::stack<Coord> ps;
    ps.push({x, y});
    do {
	  size_t x = ps.top().first, y = ps.top().second; // auto [x, y] = ps.top();
      ps.pop();
      if (grid[y*n+x])
        continue;
      grid[y*n+x] = true;
      if (x > 0)     ps.push({x - 1, y});
      if (x < n - 1) ps.push({x + 1, y});
      if (y > 0)     ps.push({x,     y - 1});
      if (y < n - 1) ps.push({x,     y + 1});
    } while (!ps.empty());
  }

  // Generates a mesh using a discretization of the domain (using a bitmap of size 2^size).
  // Assumes that domain is in [-1,1]x[-1,1].
  TriMesh regularMesh(const Point2DVector &domain, size_t size) {
    size_t n = (size_t)std::pow(2, size);
    std::vector<bool> grid(n * n, false);
    Point2D offset(-1.05, -1.05);
    double scaling = n / 2.1;

    // Init
    Point2D p0 = (domain.back() - offset) * scaling;
    int x0 = (int)p0[0], y0 = (int)p0[1];
    for (const auto &p : domain) {
      Point2D p1 = (p - offset) * scaling;
      int x1 = (int)p1[0], y1 = (int)p1[1];
      int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
      int dy = std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
      int err = (dx > dy ? dx : -dy) / 2, e2;
      while (true) {
        grid[y0*n+x0] = true;
        if (x0 == x1 && y0 == y1)
          break;
        e2 = err;
        if (e2 > -dx) { err -= dy; x0 += sx; }
        if (e2 <  dy) { err += dx; y0 += sy; }
      }
    }
    floodFill(grid, n, 0, 0);

    // Build mesh
    std::vector<size_t> row(n);
    PointVector pv;
    TriMesh result;
    size_t index = 0;

    for (size_t j = 1; j < n; ++j) {
      for (size_t i = 0; i < n - 1; ++i) {
        if (grid[j*n+i])
          continue;

        if (grid[(j-1)*n+i+1]) {
          // no NE
          if (!grid[(j-1)*n+i] && !grid[j*n+i+1])
            result.addTriangle(index, row[i], index + 1);   // N & E
        } else {
          // NE
          if (!grid[(j-1)*n+i])
            result.addTriangle(index, row[i], row[i+1]);    // N & NE
          if (!grid[j*n+i+1])
            result.addTriangle(index, row[i+1], index + 1); // E & NE
        }

        pv.emplace_back((double)i / scaling + offset[0], (double)j / scaling + offset[1], 0.0);
        row[i] = index++;
      }
    }

    result.setPoints(pv);
    return result;
  }

}

// Generates a new mesh cache using a bitmap
void
ConcaveGB::generateRegularMesh(size_t downsampling) const {
  mesh_cache_ = regularMesh(domain_, param_levels_ - downsampling);
  param_cache_.clear();
  param_cache_.reserve(mesh_cache_.points().size());
  for (const auto &p : mesh_cache_.points())
    param_cache_.emplace_back(localCoordinates(Point2D(p[0], p[1])));
}

TriMesh
ConcaveGB::evaluate(double resolution) const {
  if (last_resolution_ != resolution) {
    if (resolution > 0)
      generateDelaunayMesh(resolution);
    else
      generateRegularMesh((size_t)-resolution);
    last_resolution_ = resolution;
  }

  // Evaluate based on the cached barycentric coordinates
  def_max = 0; def_sum = 0;
  for (size_t i = 0, ie = param_cache_.size(); i != ie; ++i)
    mesh_cache_[i] = evaluate(param_cache_[i]);

  {
    std::ofstream f("deficiency.txt");
    f << "Def. max: " << def_max << std::endl;
    f << "Def. sum: " << def_sum << std::endl;
  }

  return mesh_cache_;
}

// Returns a surface point given the barycentric coordinates bc.
Point3D
ConcaveGB::evaluate(const DoubleVector &bc) const {
  size_t n = ribbons_.size();
  Point3D result(0, 0, 0);

  // Precompute the (s,d) local coordinates
  Point2DVector sds; sds.reserve(n);
  for (size_t i = 0; i < n; ++i)
    sds.push_back(barycentricSD(bc, i, parameter_dilation_));

  DoubleVector bl_s, bl_d;
  double weight_sum = 0.0;
  for (size_t side = 0; side < n; ++side) {
    size_t l = ribbons_[side].size();
    size_t d = ribbons_[side][0].size() - 1;
    bernstein(d, sds[side][0], bl_s);
    bernstein(2 * l - 1, sds[side][1], bl_d);
    for (size_t row = 0; row < l; ++row) {
      for (size_t col = 0; col <= d; ++col) {
        double blend = bl_s[col] * bl_d[row] * weight(d, sds, side, col, row);
        result += ribbons_[side][row][col] * blend;
        weight_sum += blend;
      }
    }
  }

#ifdef CONCAVE_CONTROL_POINTS
  // Add control points to the concave corners (experimental)
  for (size_t i = 0; i < n; ++i) {
    size_t ip = (i + 1) % n;
    auto &r1 = ribbons_[i];
    auto &r2 = ribbons_[ip];
    size_t l1 = r1.size(), l2 = r2.size();
    size_t d = r1[0].size() - 1;
    auto v1 = (r1[0][d] - r1[0][d-1]).normalize();
    if (v1 * (r2[1][0] - r2[0][0]) < 0)
      continue;
    Point3D cp = (r1[1][d] * 2 - r1[1][d-1] + r2[1][0] * 2 - r2[1][1]) / 2.0;
    bernstein(2 * l1 - 1, sds[i][1], bl_s);
    bernstein(2 * l2 - 1, sds[ip][1], bl_d);
    double blend = bl_s[2] * bl_d[2];
    result += cp * blend;
    weight_sum += blend;
  }
#endif

  double central_blend = 0.0;
  switch (central_weight_) {
  case CentralWeight::ORIGINAL:
    central_blend = 1.0 - weight_sum;
    break;
  case CentralWeight::ZERO:
    break;
  case CentralWeight::NTH:
    central_blend = (1.0 - weight_sum) / (double)n;
    break;
  case CentralWeight::HARMONIC:
    central_blend = n;
    for (size_t i = 0; i < n; ++i)
      central_blend *= std::pow(sds[i][1], 2);
    break;
  }

  double def = 1.0 - (weight_sum + central_blend);
  def_sum += def;
  if (std::abs(def) > std::abs(def_max))
    def_max = def;

  result += central_cp_ * central_blend;
  result /= weight_sum + central_blend;

  return result;
}

}
