#include "bezier.hh"

namespace CGB {

void bernstein(size_t n, double u, DoubleVector &coeff) {
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

Point2D bezierEval(const Point2DVector &cp, double u) {
  DoubleVector coeff;
  size_t d = cp.size() - 1;
  bernstein(d, u, coeff);
  Point2D p(0, 0);
  for (size_t j = 0; j <= d; ++j)
    p += cp[j] * coeff[j];
  return p;
}

}
