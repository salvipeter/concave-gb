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
  auto x = A.jacobiSvd(Eigen::ComputeFullV).matrixV().col(2);

  auto normal = Vector3D(x(0), x(1), x(2)).normalize();
  return { centroid, normal };
}

namespace {

  // Generates an orthogonal (u,v) coordinate system in the plane defined by `normal`.
  void localSystem(const Vector3D &normal, Vector3D &u, Vector3D &v) {
    int maxi = 0, nexti = 1;
    double max = std::abs(normal[0]), next = std::abs(normal[1]);
    if (max < next) {
      std::swap(max, next);
      std::swap(maxi, nexti);
    }
    if (std::abs(normal[2]) > max) {
      nexti = maxi;
      maxi = 2;
    } else if (std::abs(normal[2]) > next)
      nexti = 2;

    u = Vector3D(0, 0, 0);
    u[nexti] = -normal[maxi];
    u[maxi] = normal[nexti];
    u /= u.norm();
    v = normal ^ u;
  }

}

Point2DVector projectToBestFitPlane(const PointVector &pv) {
  auto lsq = fitPlane(pv);
  
  Vector3D u, v;
  localSystem(lsq.n, u, v);

  Point2DVector projected; projected.reserve(pv.size());
  for (const auto &p : pv) {
    auto q = p - lsq.n * ((p - lsq.p) * lsq.n);
    projected.emplace_back(q * u, q * v);
  }

  return projected;
}

}
