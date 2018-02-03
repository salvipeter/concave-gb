#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#include <harmonic.h>

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
  param_levels_(9), central_weight_(CentralWeight::ZERO), central_cp_(0, 0, 0),
  last_resolution_(0.0)
{
}

ConcaveGB::~ConcaveGB() {
  for (auto p : parameters_)
    harmonic_free(p);
}


void
ConcaveGB::setParamLevels(unsigned int levels) {
  param_levels_ = levels;
}

void
ConcaveGB::setCentralWeight(ConcaveGB::CentralWeight type) {
  central_weight_ = type;
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
ConcaveGB::loadControlPoints(const std::string &filename, bool generate_domain) {
  ribbons_.clear();

  size_t n, d, l;

  std::ifstream f(filename);
  f >> n;
  Point3D p;
  f >> central_cp_[0] >> central_cp_[1] >> central_cp_[2];
  ribbons_.reserve(n);
  for (size_t side = 0; side < n; ++side) {
    f >> d >> l;
    Ribbon r; r.reserve(l);
    for (size_t row = 0; row < l; ++row) {
      PointVector one_row; one_row.reserve(d + 1);
      for (size_t col = 0; col <= d; ++col) {
        f >> p[0] >> p[1] >> p[2];
        one_row.push_back(p);
      }
      r.push_back(std::move(one_row));
    }
    ribbons_.push_back(std::move(r));
  }
  
  if (generate_domain && f.good())
    generateDomain();

  return f.good();
}

namespace {

  double
  inrange(double min, double x, double max) {
    if (x < min)
      return min;
    if (x > max)
      return max;
    return x;
  }

  double
  complement_angle(double alpha) {
    if (alpha < 0)
      return -(M_PI + alpha);
    return M_PI - alpha;
  }

}

void
ConcaveGB::generateDomain() {
  // Invalidate everything
  domain_.clear();
  for (auto p : parameters_)
    harmonic_free(p);
  parameters_.clear();
  last_resolution_ = 0.0;
  param_cache_.clear();
  mesh_cache_.clear();

  size_t n = ribbons_.size();

  // Compute lengths
  DoubleVector lengths; lengths.reserve(n);
  std::transform(ribbons_.begin(), ribbons_.end(), std::back_inserter(lengths),
                 [](const auto &r) { return BCurve(r[0]).arcLength(0.0, 1.0); });
  double length_sum = std::accumulate(lengths.begin(), lengths.end(), 0.0);

  // Compute angles
  VectorVector der;
  DoubleVector angles; angles.reserve(n);
  double angle_sum = 0.0;
  size_t concave_count = 0;
  for (size_t i = 0; i < n; ++i) {
    size_t ip = (i + 1) % n;
    auto &r1 = ribbons_[i][0];
    auto &r2 = ribbons_[ip][0];
    size_t d = r1.size() - 1;
    auto v1 = (r1[d] - r1[d-1]).normalize();
    auto v2 = (r2[1] - r2[0]).normalize();
    double alpha = std::acos(inrange(-1, v1 * v2, 1));
    if (v1 * (ribbons_[ip][1][0] - r2[0]) > 0) {
      alpha *= -1;
      ++concave_count;
    }
    angles.push_back(alpha);
    angle_sum += complement_angle(angles.back());
  }
  double angle_multiplier = (n - 2 - 2 * concave_count) * M_PI / angle_sum;
  std::transform(angles.begin(), angles.end(), angles.begin(),
                 [angle_multiplier](double x) {
                   return complement_angle(complement_angle(x) * angle_multiplier);
                 });

  // Compute open domain
  domain_.resize(n);
  double dir = 0.0;
  Vector2D prev_v(0.0, 0.0);
  for (size_t i = 0; i < n; ++i) {
    domain_[i] = prev_v + Vector2D(std::cos(dir), std::sin(dir)) * lengths[i];
    dir += angles[i];
    prev_v = domain_[i];
  }

  // Compute closed domain
  double accumulated = 0.0;
  for (size_t i = 0; i < n; ++i) {
    accumulated += lengths[i];
    domain_[i] -= domain_.back() * accumulated / length_sum;
  }

  // Rescale to [-1,-1]x[1,1]
  double minx = 0.0, miny = 0.0, maxx = 0.0, maxy = 0.0;
  for (const auto &v : domain_) {
    minx = std::min(minx, v[0]); miny = std::min(miny, v[1]);
    maxx = std::max(maxx, v[0]); maxy = std::max(maxy, v[1]);
  }
  double width = std::max(maxx - minx, maxy - miny);
  Point2D topleft(-1.0, -1.0), minp(minx, miny);
  for (auto &v : domain_)
    v = topleft + (v - minp) * 2.0 / width;

  // Setup parameterization
  DoubleVector points; points.reserve(3 * n);
  for (const auto &p : domain_) {
    points.push_back(p[0]);
    points.push_back(p[1]);
    points.push_back(0);
  }
  for (size_t i = 0; i < n; ++i) {
    points[3*i+2] = 1;
    auto map = harmonic_init(n, &points[0], param_levels_, EPSILON);
    parameters_.push_back(map);
    points[3*i+2] = 0;
  }
}

namespace {

  double
  bernstein(size_t n, size_t i, double u) {
    DoubleVector tmp(n + 1, 0.0);
    tmp[n-i] = 1.0;
    double u1 = 1.0 - u;
    for (size_t k = 1; k <= n; ++k)
      for (size_t j = n; j >= k; --j)
        tmp[j] = tmp[j-1] * u + tmp[j] * u1;
    return tmp[n];
  }

}

namespace {

  size_t closest_edge(const Point2DVector &points, const Point2D &p) {
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

}

DoubleVector
ConcaveGB::localCoordinates(const Point2D &uv) const {
  size_t n = parameters_.size();
  DoubleVector result(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    double value;
    if (!harmonic_eval(parameters_[i], const_cast<double *>(uv.data()), &value)) {
      // linear interpolation
      size_t k = closest_edge(domain_, uv), km = (k + n - 1) % n;
      double x = (uv - domain_[km]).norm() / (domain_[k] - domain_[km]).norm();
      result[km] = 1 - x;
      result[k] = x;
      return result;
    }
    result[i] = value;
  }
  return result;
}

namespace {

  Point2D
  barycentricSD(const DoubleVector &bc, size_t i) {
    size_t n = bc.size(), im = (i + n - 1) % n;
    double him = bc[im], hi = bc[i];
    double d = 1.0 - him - hi;
    double s = him + hi;
    if (std::abs(s) > EPSILON)
      s = hi / s;
    return Point2D(s, d);
  }

}

double
ConcaveGB::weight(const DoubleVector &bc, size_t i, size_t j, size_t k) const {
  // TODO: could be faster if all Bernstein polynomials were computed together
  size_t n = ribbons_.size();
  size_t d = ribbons_[i][0].size() - 1;
  size_t im = (i + n - 1) % n, ip = (i + 1) % n;
  Point2D sd = barycentricSD(bc, i);
  double w = bernstein(d, j, sd[0]) * bernstein(d, k, sd[1]);
  double mu = 1.0;
  if (k < 2 && j < 2) {
    // alpha
    Point2D sdm = barycentricSD(bc, im);
    double di2 = std::pow(sd[1], 2), dim2 = std::pow(sdm[1], 2);
    double denom = dim2 + di2;
    if (denom < EPSILON)
      mu = 0.5;
    else
      mu = dim2 / denom;
  } else if (k < 2 && j > d - 2) {
    // beta
    Point2D sdp = barycentricSD(bc, ip);
    double di2 = std::pow(sd[1], 2), dip2 = std::pow(sdp[1], 2);
    double denom = dip2 + di2;
    if (denom < EPSILON)
      mu = 0.5;
    else
      mu = dip2 / denom;
  } else if (j < k || j > d - k)
    mu = 0.0;
  else if (j == k || j == d - k)
    mu = 0.5;
  return w * mu;
}

TriMesh
ConcaveGB::evaluate(double resolution) const {
  if (last_resolution_ != resolution) {
    // Generate a new mesh cache with Shewchuk's Triangle
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
    cmd << "pa" << std::fixed << resolution << "qzQ";
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
    last_resolution_ = resolution;
  }

  for (size_t i = 0, ie = param_cache_.size(); i != ie; ++i)
    mesh_cache_[i] = evaluate(param_cache_[i]);
  return mesh_cache_;
}

Point3D
ConcaveGB::evaluate(const DoubleVector &bc) const {
  size_t n = ribbons_.size();
  Point3D result(0, 0, 0);
  double weight_sum = 0.0;
  for (size_t side = 0; side < n; ++side) {
    size_t l = ribbons_[side].size();
    size_t d = ribbons_[side][0].size() - 1;
    for (size_t row = 0; row < l; ++row) {
      for (size_t col = 0; col <= d; ++col) {
        double blend = weight(bc, side, col, row);
        result += ribbons_[side][row][col] * blend;
        weight_sum += blend;
      }
    }
  }

  switch (central_weight_) {
  case CentralWeight::ORIGINAL:
    result += central_cp_ * (1.0 - weight_sum);
    break;
  case CentralWeight::ZERO:
    result /= weight_sum;
    break;
  default:
    ;
  }

  return result;
}

}
