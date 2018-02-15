#pragma once

#include <geometry.hh>

namespace LSQPlane {

using namespace Geometry;

struct Plane {
  Point3D p;
  Vector3D n;
};

Plane fitPlane(const PointVector &pv);

Point2DVector projectToBestFitPlane(const PointVector &pv);

}
