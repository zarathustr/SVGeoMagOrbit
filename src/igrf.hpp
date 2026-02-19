#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace geomag_demo {

enum class SVMode {
  None,
  Linear,
  Drift,
  LinearDrift
};

SVMode svModeFromString(const std::string& s);
std::string svModeToString(SVMode m);

class IGRFModel {
 public:
  bool load(const std::string& g_file, const std::string& h_file);

  int maxDegree() const { return nmax_; }

  // Evaluate B in spherical ECEF components [nT]: (Br, Btheta, Bphi)
  // Inputs:
  // - r_ecef_km: ECEF position in km
  // - days_since_2015: days since 2015-01-01 00:00:00 UTC
  // - N: max degree/order
  // - mode: SV handling mode
  // - drift_rate_deg_per_yr: east-positive (default: -0.2 deg/yr ~ westward drift)
  Eigen::Vector3d fieldECEF(const Eigen::Vector3d& r_ecef_km,
                            double days_since_2015,
                            int N,
                            SVMode mode,
                            double drift_rate_deg_per_yr = -0.2) const;

 private:
  int nmax_ = 0;

  // Coeff arrays indexed by [n][m] with 0<=m<=n<=nmax.
  std::vector<std::vector<double>> g0_;
  std::vector<std::vector<double>> gsv_;
  std::vector<std::vector<double>> h0_;
  std::vector<std::vector<double>> hsv_;

  double coeffG(int n, int m, double days_since_2015, bool use_linear_sv) const;
  double coeffH(int n, int m, double days_since_2015, bool use_linear_sv) const;
};

}  // namespace geomag_demo
