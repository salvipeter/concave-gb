#pragma once

#include <geometry.hh>

namespace CGB {

using namespace Geometry;

void bernstein(size_t n, double u, DoubleVector &coeff);

Point2D bezierEval(const Point2DVector &cp, double u);

}
