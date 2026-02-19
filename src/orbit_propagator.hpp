#pragma once

#include <Eigen/Dense>

namespace geomag_demo {

struct OrbitParams {
  double mu_km3_s2 = 398600.4415;        // Earth's GM [km^3/s^2]
  double Re_km = 6378.1363;             // Earth radius [km]
  double J2 = 1.08262668e-3;            // J2
};

// One-step RK4 propagation of a 6D ECI/J2000 state [m, m/s], with 2-body + J2.
// If apply_axis_swap is true, reproduces the legacy STK axis rotation used in the MATLAB demo.
Eigen::Matrix<double, 6, 1> propagateRK4(const Eigen::Matrix<double, 6, 1>& x_m,
                                        double dt_s,
                                        const OrbitParams& p,
                                        bool apply_axis_swap = true);

}  // namespace geomag_demo
