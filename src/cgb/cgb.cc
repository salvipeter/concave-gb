#include <cassert>
#include <fstream>

#include "cgb.hh"

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
#include <triangle.h>

namespace CGB {

ConcaveGB::ConcaveGB() : n_(0), d_(0) {
}

ConcaveGB::~ConcaveGB() {
  for (auto p : parameters_)
    harmonic_free(p);
}

void
ConcaveGB::loadControlPoints(const std::string &filename) {
}

void
ConcaveGB::loadControlPoints(size_t n, size_t d, std::vector<Point3D>) {
}

void
ConcaveGB::generateDomain() {
  assert(n_ > 0);
}

TriMesh
ConcaveGB::evaluate(double resolution) const {
  return TriMesh();
}

Point3D
ConcaveGB::evaluate(const Point2D &uv) const {
  return Point3D(0,0,0);
}

}
