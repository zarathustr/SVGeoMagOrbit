#include <Eigen/Dense>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "eci_ecef.hpp"
#include "igrf.hpp"
#include "orbit_propagator.hpp"
#include "stk_loader.hpp"
#include "time_utils.hpp"

#ifndef GEOMAG_DEMO_DATA_DIR
#define GEOMAG_DEMO_DATA_DIR ""
#endif

namespace fs = std::filesystem;
using geomag_demo::DateTime;
using geomag_demo::IGRFModel;
using geomag_demo::OrbitParams;
using geomag_demo::SVMode;
using geomag_demo::StkRow;

struct Args {
  fs::path data_dir;
  std::string posvel_file = "Satellite2 J2000 Position Velocity.txt";

  std::string g_file = "igrfSg2015.txt";
  std::string h_file = "igrfSh2015.txt";

  std::size_t steps = 1200;         // default: 20 minutes at 1 Hz
  double dt_s = 1.0;

  int igrf_degree_full = 12;
  int igrf_degree_dipole = 1;

  SVMode truth_sv_mode = SVMode::Linear;
  SVMode est_sv_mode = SVMode::Linear;
  double drift_rate_deg_per_yr = -0.2;  // east-positive

  double sigma_B_nT = 10.0;

  // Initial covariance / perturbations
  double init_pos_sigma_m = 100e3;   // 100 km
  double init_vel_sigma_mps = 10.0;  // 10 m/s

  bool apply_axis_swap = true;
  unsigned int seed = 1;

  bool write_csv = true;
  std::string csv_file = "geomag_demo_results.csv";
};

static void printUsage(const char* exe) {
  std::cout
      << "Usage: " << exe << " [options]\n\n"
      << "Options:\n"
      << "  --data_dir <path>        Directory containing STK Satellite2*.txt files\n"
      << "  --posvel_file <name>     Position/velocity report filename (default: \"Satellite2 J2000 Position Velocity.txt\")\n"
      << "  --steps <N>              Number of samples to process (default: 1200)\n"
      << "  --dt <sec>               Sample interval in seconds (default: 1)\n"
      << "  --sv_mode <mode>         Estimator SV mode: none|linear|drift|linear_drift (default: linear)\n"
      << "  --truth_sv_mode <mode>   Truth SV mode used to generate measurements (default: linear)\n"
      << "  --drift_rate <deg/yr>    Drift rate (east-positive; default: -0.2)\n"
      << "  --sigma_B <nT>           Magnetometer 1-sigma noise per axis (default: 10)\n"
      << "  --init_pos_sigma <m>     Initial position 1-sigma perturbation (default: 100000)\n"
      << "  --init_vel_sigma <m/s>   Initial velocity 1-sigma perturbation (default: 10)\n"
      << "  --no_axis_swap           Disable the legacy STK axis rotation used in the MATLAB demo\n"
      << "  --no_csv                 Do not write CSV output\n"
      << "  --csv <file>             CSV output filename (default: geomag_demo_results.csv)\n"
      << "  --seed <int>             RNG seed (default: 1)\n"
      << "  -h, --help               Show this help\n";
}

static bool parseArgs(int argc, char** argv, Args& a) {
  // Default data dir from compile definition (if provided)
  if (std::string(GEOMAG_DEMO_DATA_DIR).size() > 0) {
    a.data_dir = fs::path(GEOMAG_DEMO_DATA_DIR);
  } else {
    a.data_dir = fs::current_path();
  }

  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    auto needValue = [&](const std::string& k) -> const char* {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << k << "\n";
        return nullptr;
      }
      return argv[++i];
    };

    if (key == "-h" || key == "--help") {
      return false;
    } else if (key == "--data_dir") {
      const char* v = needValue(key);
      if (!v) return false;
      a.data_dir = fs::path(v);
    } else if (key == "--posvel_file") {
      const char* v = needValue(key);
      if (!v) return false;
      a.posvel_file = v;
    } else if (key == "--steps") {
      const char* v = needValue(key);
      if (!v) return false;
      a.steps = static_cast<std::size_t>(std::stoll(v));
    } else if (key == "--dt") {
      const char* v = needValue(key);
      if (!v) return false;
      a.dt_s = std::stod(v);
    } else if (key == "--sv_mode") {
      const char* v = needValue(key);
      if (!v) return false;
      a.est_sv_mode = geomag_demo::svModeFromString(v);
    } else if (key == "--truth_sv_mode") {
      const char* v = needValue(key);
      if (!v) return false;
      a.truth_sv_mode = geomag_demo::svModeFromString(v);
    } else if (key == "--drift_rate") {
      const char* v = needValue(key);
      if (!v) return false;
      a.drift_rate_deg_per_yr = std::stod(v);
    } else if (key == "--sigma_B") {
      const char* v = needValue(key);
      if (!v) return false;
      a.sigma_B_nT = std::stod(v);
    } else if (key == "--init_pos_sigma") {
      const char* v = needValue(key);
      if (!v) return false;
      a.init_pos_sigma_m = std::stod(v);
    } else if (key == "--init_vel_sigma") {
      const char* v = needValue(key);
      if (!v) return false;
      a.init_vel_sigma_mps = std::stod(v);
    } else if (key == "--no_axis_swap") {
      a.apply_axis_swap = false;
    } else if (key == "--no_csv") {
      a.write_csv = false;
    } else if (key == "--csv") {
      const char* v = needValue(key);
      if (!v) return false;
      a.csv_file = v;
      a.write_csv = true;
    } else if (key == "--seed") {
      const char* v = needValue(key);
      if (!v) return false;
      a.seed = static_cast<unsigned int>(std::stoul(v));
    } else {
      std::cerr << "Unknown option: " << key << "\n";
      return false;
    }
  }

  return true;
}

static Eigen::Vector3d magneticFieldECEF_sph(const IGRFModel& igrf,
                                            const Eigen::Matrix<double, 6, 1>& x_eci_m,
                                            double gast_deg,
                                            double days_since_2015,
                                            int degree,
                                            SVMode mode,
                                            double drift_rate_deg_per_yr) {
  const Eigen::Vector3d r_eci_m = x_eci_m.segment<3>(0);
  const Eigen::Vector3d r_ecef_m = geomag_demo::eciToEcefPosition(r_eci_m, gast_deg);
  const Eigen::Vector3d r_ecef_km = r_ecef_m * 1e-3;
  return igrf.fieldECEF(r_ecef_km, days_since_2015, degree, mode, drift_rate_deg_per_yr);
}

static Eigen::Matrix<double, 6, 6> numericalJacobianF(
    const Eigen::Matrix<double, 6, 1>& x,
    double dt,
    const OrbitParams& p,
    bool apply_axis_swap) {
  Eigen::Matrix<double, 6, 6> F;
  F.setZero();

  const double delta_pos = 10.0;   // [m]
  const double delta_vel = 0.01;   // [m/s]

  for (int i = 0; i < 6; ++i) {
    const double d = (i < 3) ? delta_pos : delta_vel;
    Eigen::Matrix<double, 6, 1> xp = x;
    Eigen::Matrix<double, 6, 1> xm = x;
    xp(i) += d;
    xm(i) -= d;

    const Eigen::Matrix<double, 6, 1> fp = geomag_demo::propagateRK4(xp, dt, p, apply_axis_swap);
    const Eigen::Matrix<double, 6, 1> fm = geomag_demo::propagateRK4(xm, dt, p, apply_axis_swap);

    F.col(i) = (fp - fm) / (2.0 * d);
  }

  return F;
}

static Eigen::Matrix<double, 3, 6> numericalJacobianH_dipole(
    const IGRFModel& igrf,
    const Eigen::Matrix<double, 6, 1>& x,
    double gast_deg,
    double days_since_2015,
    SVMode mode,
    double drift_rate_deg_per_yr) {
  Eigen::Matrix<double, 3, 6> H;
  H.setZero();

  const double delta_pos = 10.0;  // [m]

  for (int i = 0; i < 3; ++i) {
    Eigen::Matrix<double, 6, 1> xp = x;
    Eigen::Matrix<double, 6, 1> xm = x;
    xp(i) += delta_pos;
    xm(i) -= delta_pos;

    const Eigen::Vector3d hp = magneticFieldECEF_sph(igrf, xp, gast_deg, days_since_2015, 1, mode, drift_rate_deg_per_yr);
    const Eigen::Vector3d hm = magneticFieldECEF_sph(igrf, xm, gast_deg, days_since_2015, 1, mode, drift_rate_deg_per_yr);

    H.col(i) = (hp - hm) / (2.0 * delta_pos);
  }

  return H;
}

int main(int argc, char** argv) {
  Args args;
  if (!parseArgs(argc, argv, args)) {
    printUsage(argv[0]);
    return (argc > 1) ? 1 : 0;
  }

  const fs::path posvel_path = args.data_dir / args.posvel_file;

  std::vector<StkRow> rows;
  if (!geomag_demo::loadStkReport(posvel_path, rows, args.steps)) {
    std::cerr << "Failed to load STK report: " << posvel_path << "\n";
    std::cerr << "Tip: pass --data_dir /path/to/STK_exports\n";
    return 2;
  }

  const std::size_t N = std::min<std::size_t>(args.steps, rows.size());
  if (N < 2) {
    std::cerr << "Not enough samples in the report.\n";
    return 3;
  }

  // Load IGRF coefficients (expect files in the working directory, or supply absolute paths).
  IGRFModel igrf;
  if (!igrf.load(args.g_file, args.h_file)) {
    std::cerr << "Failed to load IGRF coefficient files: " << args.g_file << ", " << args.h_file << "\n";
    std::cerr << "Tip: run from the build directory where CMake copies igrfSg2015.txt / igrfSh2015.txt\n";
    return 4;
  }

  std::cout << "Loaded " << N << " samples from: " << posvel_path << "\n";
  std::cout << "Estimator SV mode: " << geomag_demo::svModeToString(args.est_sv_mode)
            << ", Truth SV mode: " << geomag_demo::svModeToString(args.truth_sv_mode)
            << ", drift_rate=" << args.drift_rate_deg_per_yr << " deg/yr\n";

  // Parse truth states
  std::vector<Eigen::Matrix<double, 6, 1>> x_true(N);
  std::vector<double> gast_deg(N);
  std::vector<double> days2015(N);

  for (std::size_t k = 0; k < N; ++k) {
    if (rows[k].values.size() < 6) {
      std::cerr << "Row " << k << " has <6 numeric columns; cannot interpret as pos/vel.\n";
      return 5;
    }

    // STK report values are typically [km, km, km, km/s, km/s, km/s]. Convert to [m, m/s].
    Eigen::Matrix<double, 6, 1> x;
    for (int i = 0; i < 6; ++i) {
      x(i) = rows[k].values[i] * 1e3;
    }
    x_true[k] = x;

    DateTime dt;
    if (!geomag_demo::parseStkTimestamp(rows[k].timestamp, dt)) {
      std::cerr << "Failed to parse timestamp at row " << k << ": '" << rows[k].timestamp << "'\n";
      return 6;
    }

    const double jd = geomag_demo::toJulianDateUTC(dt);
    gast_deg[k] = geomag_demo::jdToGAST_deg(jd);
    days2015[k] = geomag_demo::daysSince2015(jd);
  }

  // Generate synthetic magnetometer measurements (3-axis in spherical ECEF)
  std::mt19937 rng(args.seed);
  std::normal_distribution<double> n01(0.0, 1.0);

  std::vector<Eigen::Vector3d> B_meas(N);
  for (std::size_t k = 0; k < N; ++k) {
    const Eigen::Vector3d B_true = magneticFieldECEF_sph(
        igrf, x_true[k], gast_deg[k], days2015[k], args.igrf_degree_full, args.truth_sv_mode, args.drift_rate_deg_per_yr);

    Eigen::Vector3d noise;
    noise << args.sigma_B_nT * n01(rng), args.sigma_B_nT * n01(rng), args.sigma_B_nT * n01(rng);
    B_meas[k] = B_true + noise;
  }

  // EKF init
  OrbitParams orb;

  Eigen::Matrix<double, 6, 1> x_est = x_true[0];
  for (int i = 0; i < 3; ++i) x_est(i) += args.init_pos_sigma_m * n01(rng);
  for (int i = 3; i < 6; ++i) x_est(i) += args.init_vel_sigma_mps * n01(rng);

  Eigen::Matrix<double, 6, 6> P = Eigen::Matrix<double, 6, 6>::Zero();
  for (int i = 0; i < 3; ++i) P(i, i) = args.init_pos_sigma_m * args.init_pos_sigma_m;
  for (int i = 3; i < 6; ++i) P(i, i) = args.init_vel_sigma_mps * args.init_vel_sigma_mps;

  // Process noise (demo-level tuning)
  Eigen::Matrix<double, 6, 6> Q = Eigen::Matrix<double, 6, 6>::Zero();
  const double q_pos = 1.0;   // [m^2]
  const double q_vel = 1e-2;  // [(m/s)^2]
  for (int i = 0; i < 3; ++i) Q(i, i) = q_pos;
  for (int i = 3; i < 6; ++i) Q(i, i) = q_vel;

  // Measurement noise
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity() * (args.sigma_B_nT * args.sigma_B_nT);

  // Logs
  std::vector<double> pos_err_km(N);
  std::vector<double> vel_err_ms(N);

  // Optional CSV output
  std::ofstream csv;
  if (args.write_csv) {
    csv.open(args.csv_file);
    if (!csv.is_open()) {
      std::cerr << "Warning: could not open CSV file for writing: " << args.csv_file << "\n";
    } else {
      csv << "k,time_days_since2015,pos_err_km,vel_err_ms,";
      csv << "x_true_m,y_true_m,z_true_m,vx_true_mps,vy_true_mps,vz_true_mps,";
      csv << "x_est_m,y_est_m,z_est_m,vx_est_mps,vy_est_mps,vz_est_mps\n";
    }
  }

  // EKF loop
  for (std::size_t k = 0; k < N; ++k) {
    Eigen::Matrix<double, 6, 1> x_pred = x_est;
    Eigen::Matrix<double, 6, 6> P_pred = P;

    if (k > 0) {
      x_pred = geomag_demo::propagateRK4(x_est, args.dt_s, orb, args.apply_axis_swap);
      const Eigen::Matrix<double, 6, 6> F = numericalJacobianF(x_est, args.dt_s, orb, args.apply_axis_swap);
      P_pred = F * P * F.transpose() + Q;
    }

    // Predicted full and dipole fields at the predicted state
    const Eigen::Vector3d B_full_pred = magneticFieldECEF_sph(
        igrf, x_pred, gast_deg[k], days2015[k], args.igrf_degree_full, args.est_sv_mode, args.drift_rate_deg_per_yr);

    const Eigen::Vector3d B_dip_pred = magneticFieldECEF_sph(
        igrf, x_pred, gast_deg[k], days2015[k], args.igrf_degree_dipole, args.est_sv_mode, args.drift_rate_deg_per_yr);

    const Eigen::Vector3d B_hot_pred = B_full_pred - B_dip_pred;

    // The dipole-optimized measurement equation:
    //   y = B_meas - B_hot_pred
    //   h(x) = B_dip(x)
    // innovation = y - h(x_pred)
    const Eigen::Vector3d y = B_meas[k] - B_hot_pred;
    const Eigen::Vector3d innov = y - B_dip_pred;

    const Eigen::Matrix<double, 3, 6> H = numericalJacobianH_dipole(
        igrf, x_pred, gast_deg[k], days2015[k], args.est_sv_mode, args.drift_rate_deg_per_yr);

    const Eigen::Matrix3d S = H * P_pred * H.transpose() + R;
    const Eigen::Matrix<double, 6, 3> K = P_pred * H.transpose() * S.inverse();

    x_est = x_pred + K * innov;

    const Eigen::Matrix<double, 6, 6> I = Eigen::Matrix<double, 6, 6>::Identity();
    const Eigen::Matrix<double, 6, 6> KH = K * H;
    // Joseph form for numerical stability
    P = (I - KH) * P_pred * (I - KH).transpose() + K * R * K.transpose();

    // Errors vs truth
    const Eigen::Vector3d e_pos = x_true[k].segment<3>(0) - x_est.segment<3>(0);
    const Eigen::Vector3d e_vel = x_true[k].segment<3>(3) - x_est.segment<3>(3);

    pos_err_km[k] = e_pos.norm() * 1e-3;
    vel_err_ms[k] = e_vel.norm();

    if (csv.is_open()) {
      csv << k << "," << days2015[k] << "," << pos_err_km[k] << "," << vel_err_ms[k] << ",";
      for (int i = 0; i < 6; ++i) csv << x_true[k](i) << (i == 5 ? "," : ",");
      for (int i = 0; i < 6; ++i) {
        csv << x_est(i);
        if (i < 5) csv << ",";
      }
      csv << "\n";
    }
  }

  if (csv.is_open()) {
    csv.close();
    std::cout << "Wrote CSV: " << args.csv_file << "\n";
  }

  // Compute RMS
  double pos_rms2 = 0.0;
  double vel_rms2 = 0.0;
  for (std::size_t k = 0; k < N; ++k) {
    pos_rms2 += pos_err_km[k] * pos_err_km[k];
    vel_rms2 += vel_err_ms[k] * vel_err_ms[k];
  }
  pos_rms2 /= static_cast<double>(N);
  vel_rms2 /= static_cast<double>(N);

  std::cout << "\n==== Demo summary ====\n";
  std::cout << "Samples: " << N << ", dt=" << args.dt_s << " s\n";
  std::cout << "Position RMS: " << std::sqrt(pos_rms2) << " km\n";
  std::cout << "Velocity RMS: " << std::sqrt(vel_rms2) << " m/s\n";
  std::cout << "======================\n";

  return 0;
}
