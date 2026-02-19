#include "orbit_propagator.hpp"

#include <cmath>

namespace geomag_demo {

static Eigen::Matrix<double, 6, 1> eomTwoBodyJ2(const Eigen::Matrix<double, 6, 1>& y_km,
                                                const OrbitParams& p) {
  Eigen::Vector3d r = y_km.segment<3>(0);
  Eigen::Vector3d v = y_km.segment<3>(3);

  const double r2 = r.squaredNorm();
  const double r1 = std::sqrt(r2);
  const double r3 = r2 * r1;

  Eigen::Vector3d a = -p.mu_km3_s2 * r / r3;

  // J2
  const double z2 = r.z() * r.z();
  const double r5 = r3 * r2;
  const double factor = 1.5 * p.J2 * p.mu_km3_s2 * (p.Re_km * p.Re_km) / r5;
  const double f = 5.0 * z2 / r2;

  Eigen::Vector3d aJ2;
  aJ2 << r.x() * (f - 1.0),
         r.y() * (f - 1.0),
         r.z() * (f - 3.0);

  a += factor * aJ2;

  Eigen::Matrix<double, 6, 1> ydot;
  ydot.segment<3>(0) = v;
  ydot.segment<3>(3) = a;
  return ydot;
}

static Eigen::Matrix<double, 6, 1> applyAxisSwap(const Eigen::Matrix<double, 6, 1>& x_km) {
  // Equivalent to MATLAB:
  // y = [ y(2); -y(1); y(3); y(5); -y(4); y(6) ];
  Eigen::Matrix<double, 6, 1> y;
  y << x_km(1),
       -x_km(0),
       x_km(2),
       x_km(4),
       -x_km(3),
       x_km(5);
  return y;
}

static Eigen::Matrix<double, 6, 1> inverseAxisSwap(const Eigen::Matrix<double, 6, 1>& y_km) {
  // Equivalent to MATLAB:
  // y_next = [ -y_next(2);
  //             y_next(1);
  //             y_next(3);
  //           - y_next(5);
  //             y_next(4);
  //             y_next(6) ];
  Eigen::Matrix<double, 6, 1> x;
  x << -y_km(1),
        y_km(0),
        y_km(2),
       -y_km(4),
        y_km(3),
        y_km(5);
  return x;
}

Eigen::Matrix<double, 6, 1> propagateRK4(const Eigen::Matrix<double, 6, 1>& x_m,
                                        double dt_s,
                                        const OrbitParams& p,
                                        bool apply_axis_swap) {
  // meters -> km
  Eigen::Matrix<double, 6, 1> y = x_m * 1e-3;

  if (apply_axis_swap) {
    y = applyAxisSwap(y);
  }

  const Eigen::Matrix<double, 6, 1> k1 = eomTwoBodyJ2(y, p);
  const Eigen::Matrix<double, 6, 1> k2 = eomTwoBodyJ2(y + 0.5 * dt_s * k1, p);
  const Eigen::Matrix<double, 6, 1> k3 = eomTwoBodyJ2(y + 0.5 * dt_s * k2, p);
  const Eigen::Matrix<double, 6, 1> k4 = eomTwoBodyJ2(y + dt_s * k3, p);

  Eigen::Matrix<double, 6, 1> y_next = y + (dt_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

  if (apply_axis_swap) {
    y_next = inverseAxisSwap(y_next);
  }

  // km -> m
  return y_next * 1e3;
}

}  // namespace geomag_demo
