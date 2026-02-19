#include "eci_ecef.hpp"

#include <cmath>

namespace geomag_demo {

static constexpr double kPi = 3.14159265358979323846;
static double deg2rad(double deg) { return deg * kPi / 180.0; }

Eigen::Matrix3d R3_gast(double gast_deg) {
  const double th = deg2rad(gast_deg);
  const double c = std::cos(th);
  const double s = std::sin(th);

  Eigen::Matrix3d R;
  R <<  c,  s, 0,
       -s,  c, 0,
        0,  0, 1;
  return R;
}

Eigen::Vector3d eciToEcefPosition(const Eigen::Vector3d& r_eci_m, double gast_deg) {
  return R3_gast(gast_deg) * r_eci_m;
}

}  // namespace geomag_demo
