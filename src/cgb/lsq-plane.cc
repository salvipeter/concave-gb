#include "Eigen/SVD"

#include "lsq-plane.hh"

namespace LSQPlane {

Plane fitPlane(const PointVector &pv) {
  // Compute centroid
  size_t n = pv.size();
  Point3D centroid(0, 0, 0);
  for (const auto &p : pv)
    centroid += p;
  centroid /= n;

  // Solve by singular value decomposition
  Eigen::MatrixXd A(n, 3);
  for (size_t i = 0; i < n; ++i) {
    auto p = pv[i] - centroid;
    A.row(i) << p[0], p[1], p[2];
  }
  auto svd = A.jacobiSvd(Eigen::ComputeThinV);

  // Extract principal directions
  auto vec = [](const auto &v) { return Vector3D(v(0), v(1), v(2)); };
  auto u = vec(svd.matrixV().col(0));
  auto v = vec(svd.matrixV().col(1));
  auto normal = vec(svd.matrixV().col(2));

  return { centroid, u, v, normal };
}

Point2DVector projectToBestFitPlane(const PointVector &pv) {
  auto lsq = fitPlane(pv);
  
  Point2DVector projected; projected.reserve(pv.size());
  for (const auto &p : pv) {
    auto q = p - lsq.n * ((p - lsq.p) * lsq.n);
    projected.emplace_back(q * lsq.u, q * lsq.v);
  }

  return projected;
}

}
