#pragma once

#include <Eigen/Dense>

namespace geomag_demo {

// ECI(J2000) -> ECEF rotation matrix using GAST angle (deg).
Eigen::Matrix3d R3_gast(double gast_deg);

// Convert ECI position [m] to ECEF position [m] using GAST.
Eigen::Vector3d eciToEcefPosition(const Eigen::Vector3d& r_eci_m, double gast_deg);

}  // namespace geomag_demo
