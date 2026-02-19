#include "igrf.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <sstream>
#include <tuple>

namespace geomag_demo {

static constexpr double kPi = 3.14159265358979323846;
static double deg2rad(double deg) { return deg * kPi / 180.0; }

SVMode svModeFromString(const std::string& s) {
  std::string t;
  t.reserve(s.size());
  for (char c : s) {
    t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  }
  if (t == "none") return SVMode::None;
  if (t == "linear") return SVMode::Linear;
  if (t == "drift") return SVMode::Drift;
  if (t == "linear_drift" || t == "linear+drift" || t == "linear-drift") return SVMode::LinearDrift;
  return SVMode::Linear;  // default
}

std::string svModeToString(SVMode m) {
  switch (m) {
    case SVMode::None: return "none";
    case SVMode::Linear: return "linear";
    case SVMode::Drift: return "drift";
    case SVMode::LinearDrift: return "linear_drift";
    default: return "linear";
  }
}

static bool readCoeffFile(const std::string& path,
                          std::vector<std::tuple<int,int,double,double>>& out) {
  out.clear();
  std::ifstream ifs(path);
  if (!ifs.is_open()) return false;

  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty()) continue;
    // allow tab or space separated
    std::istringstream iss(line);
    int n = 0, m = 0;
    double v = 0.0, sv = 0.0;
    if (!(iss >> n >> m >> v >> sv)) {
      continue;
    }
    out.emplace_back(n, m, v, sv);
  }
  return !out.empty();
}

bool IGRFModel::load(const std::string& g_file, const std::string& h_file) {
  std::vector<std::tuple<int,int,double,double>> g_rows;
  std::vector<std::tuple<int,int,double,double>> h_rows;
  if (!readCoeffFile(g_file, g_rows)) return false;
  if (!readCoeffFile(h_file, h_rows)) return false;

  int nmax = 0;
  for (const auto& r : g_rows) nmax = std::max(nmax, std::get<0>(r));
  for (const auto& r : h_rows) nmax = std::max(nmax, std::get<0>(r));
  if (nmax <= 0) return false;

  nmax_ = nmax;
  g0_.assign(nmax_ + 1, std::vector<double>(nmax_ + 1, 0.0));
  gsv_.assign(nmax_ + 1, std::vector<double>(nmax_ + 1, 0.0));
  h0_.assign(nmax_ + 1, std::vector<double>(nmax_ + 1, 0.0));
  hsv_.assign(nmax_ + 1, std::vector<double>(nmax_ + 1, 0.0));

  for (const auto& r : g_rows) {
    const int n = std::get<0>(r);
    const int m = std::get<1>(r);
    if (n < 0 || n > nmax_ || m < 0 || m > nmax_) continue;
    g0_[n][m] = std::get<2>(r);
    gsv_[n][m] = std::get<3>(r);
  }
  for (const auto& r : h_rows) {
    const int n = std::get<0>(r);
    const int m = std::get<1>(r);
    if (n < 0 || n > nmax_ || m < 0 || m > nmax_) continue;
    h0_[n][m] = std::get<2>(r);
    hsv_[n][m] = std::get<3>(r);
  }

  return true;
}

double IGRFModel::coeffG(int n, int m, double days_since_2015, bool use_linear_sv) const {
  const double years = days_since_2015 / 365.0;
  return g0_[n][m] + (use_linear_sv ? (gsv_[n][m] * years) : 0.0);
}

double IGRFModel::coeffH(int n, int m, double days_since_2015, bool use_linear_sv) const {
  const double years = days_since_2015 / 365.0;
  return h0_[n][m] + (use_linear_sv ? (hsv_[n][m] * years) : 0.0);
}

Eigen::Vector3d IGRFModel::fieldECEF(const Eigen::Vector3d& r_ecef_km,
                                    double days_since_2015,
                                    int N,
                                    SVMode mode,
                                    double drift_rate_deg_per_yr) const {
  const double x = r_ecef_km.x();
  const double y = r_ecef_km.y();
  const double z = r_ecef_km.z();

  const double r2 = x * x + y * y + z * z;
  if (r2 <= 0.0) {
    return Eigen::Vector3d::Zero();
  }

  const double r = std::sqrt(r2);
  const double phi0 = std::atan2(y, x);          // longitude [rad]
  const double theta = std::acos(z / r);         // colatitude [rad]

  const bool use_linear_sv = (mode == SVMode::Linear || mode == SVMode::LinearDrift);
  const bool use_drift = (mode == SVMode::Drift || mode == SVMode::LinearDrift);

  double phi = phi0;
  if (use_drift) {
    const double years = days_since_2015 / 365.0;
    phi = phi0 + deg2rad(drift_rate_deg_per_yr * years);
  }

  int Nuse = std::min(N, nmax_);
  if (Nuse < 1) Nuse = 1;

  // Reference radius used by IGRF [km]
  constexpr double a = 6371.2;

  const double sin_t = std::sin(theta);
  const double cos_t = std::cos(theta);
  const double sin_t_safe = (std::abs(sin_t) < 1e-12) ? (sin_t >= 0 ? 1e-12 : -1e-12) : sin_t;

  double Br = 0.0;
  double Bt = 0.0;
  double Bp = 0.0;

  // Recursive associated Legendre (matches the compact MATLAB implementation)
  double P11 = 1.0;
  double P10 = P11;
  double P20 = 0.0;

  double dP11 = 0.0;
  double dP10 = dP11;
  double dP20 = 0.0;

  for (int m = 0; m <= Nuse; ++m) {
    for (int n = 1; n <= Nuse; ++n) {
      if (m > n) continue;

      double P2 = 0.0;
      double dP2 = 0.0;

      if (n == m) {
        P2 = sin_t * P11;
        dP2 = sin_t * dP11 + cos_t * P11;

        P11 = P2;
        P10 = P11;
        P20 = 0.0;

        dP11 = dP2;
        dP10 = dP11;
        dP20 = 0.0;
      } else if (n == 1) {
        P2 = cos_t * P10;
        dP2 = cos_t * dP10 - sin_t * P10;

        P20 = P10;
        P10 = P2;

        dP20 = dP10;
        dP10 = dP2;
      } else {
        const double K = ((n - 1.0) * (n - 1.0) - m * m) / ((2.0 * n - 1.0) * (2.0 * n - 3.0));
        P2 = cos_t * P10 - K * P20;
        dP2 = cos_t * dP10 - sin_t * P10 - K * dP20;

        P20 = P10;
        P10 = P2;

        dP20 = dP10;
        dP10 = dP2;
      }

      const double gnm = coeffG(n, m, days_since_2015, use_linear_sv);
      const double hnm = coeffH(n, m, days_since_2015, use_linear_sv);

      const double cos_mphi = std::cos(m * phi);
      const double sin_mphi = std::sin(m * phi);
      const double term = gnm * cos_mphi + hnm * sin_mphi;

      const double ar = a / r;
      const double ar_n2 = std::pow(ar, n + 2);

      Br += ar_n2 * (n + 1.0) * term * P2;
      Bt -= ar_n2 * term * dP2;

      if (m != 0) {
        const double term_phi = m * (-gnm * sin_mphi + hnm * cos_mphi) * P2;
        Bp -= ar_n2 * term_phi / sin_t_safe;
      }
    }
  }

  return Eigen::Vector3d(Br, Bt, Bp);
}

}  // namespace geomag_demo
