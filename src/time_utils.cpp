#include "time_utils.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <regex>
#include <unordered_map>

namespace geomag_demo {

static constexpr double kPi = 3.14159265358979323846;

static std::string trim(const std::string& s) {
  auto is_ws = [](unsigned char c) { return std::isspace(c) != 0; };
  auto b = std::find_if_not(s.begin(), s.end(), is_ws);
  auto e = std::find_if_not(s.rbegin(), s.rend(), is_ws).base();
  if (b >= e) return std::string();
  return std::string(b, e);
}

static int monthFromAbbrev(const std::string& m) {
  static const std::unordered_map<std::string, int> kMap = {
      {"JAN", 1}, {"FEB", 2}, {"MAR", 3}, {"APR", 4}, {"MAY", 5}, {"JUN", 6},
      {"JUL", 7}, {"AUG", 8}, {"SEP", 9}, {"OCT", 10}, {"NOV", 11}, {"DEC", 12}};

  std::string u;
  u.reserve(m.size());
  for (char c : m) u.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));

  auto it = kMap.find(u);
  return (it == kMap.end()) ? 0 : it->second;
}

bool parseStkTimestamp(const std::string& ts, DateTime& out) {
  const std::string s = trim(ts);

  // Typical STK UTCG timestamp is fixed-width 24 chars, e.g.
  // "01 Mar 2019 00:00:00.000"
  // Some STK exports may omit leading zeros for the hour.
  static const std::regex re(
      R"(^\s*(\d{1,2})\s+([A-Za-z]{3})\s+(\d{4})\s+(\d{1,2}):(\d{2}):(\d{2}(?:\.\d+)?)\s*$)");

  std::smatch m;
  if (!std::regex_match(s, m, re)) {
    return false;
  }

  const int day = std::stoi(m[1].str());
  const int month = monthFromAbbrev(m[2].str());
  const int year = std::stoi(m[3].str());
  const int hour = std::stoi(m[4].str());
  const int minute = std::stoi(m[5].str());
  const double second = std::stod(m[6].str());

  if (month < 1 || month > 12) return false;
  if (day < 1 || day > 31) return false;

  out.year = year;
  out.month = month;
  out.day = day;
  out.hour = hour;
  out.minute = minute;
  out.second = second;
  return true;
}

// Gregorian calendar to Julian date (UTC).
// Reference: common astronomical conversion.
double toJulianDateUTC(const DateTime& t) {
  int Y = t.year;
  int M = t.month;
  const int D = t.day;

  if (M <= 2) {
    Y -= 1;
    M += 12;
  }

  const int A = static_cast<int>(std::floor(Y / 100.0));
  const int B = 2 - A + static_cast<int>(std::floor(A / 4.0));

  const double day_fraction = (t.hour + (t.minute + t.second / 60.0) / 60.0) / 24.0;

  const double JD = std::floor(365.25 * (Y + 4716)) + std::floor(30.6001 * (M + 1)) + D + B - 1524.5 + day_fraction;
  return JD;
}

static double wrap360(double deg) {
  double x = std::fmod(deg, 360.0);
  if (x < 0) x += 360.0;
  return x;
}

static double sind(double deg) { return std::sin(deg * kPi / 180.0); }
static double cosd(double deg) { return std::cos(deg * kPi / 180.0); }

// MATLAB JD2GMST equivalent.
double jdToGMST_deg(double jd) {
  // Julian date of the previous midnight (JD0)
  const double JDmin = std::floor(jd) - 0.5;
  const double JDmax = std::floor(jd) + 0.5;
  double JD0 = JDmin;
  if (jd > JDmax) {
    JD0 = JDmax;
  }

  const double H = (jd - JD0) * 24.0;          // hours past previous midnight
  const double D = jd - 2451545.0;            // days since J2000
  const double D0 = JD0 - 2451545.0;
  const double T = D / 36525.0;               // centuries since J2000

  const double GMST_hours = std::fmod(6.697374558 + 0.06570982441908 * D0 + 1.00273790935 * H + 0.000026 * (T * T), 24.0);
  double gmst = GMST_hours * 15.0;            // degrees
  gmst = wrap360(gmst);
  return gmst;
}

// MATLAB JD2GAST equivalent.
double jdToGAST_deg(double jd) {
  const double THETAm = jdToGMST_deg(jd);
  const double T = (jd - 2451545.0) / 36525.0;

  const double EPSILONm = 23.439291 - 0.0130111 * T - 1.64e-7 * (T * T) + 5.04e-7 * (T * T * T);

  const double L = 280.4665 + 36000.7698 * T;
  const double dL = 218.3165 + 481267.8813 * T;
  const double OMEGA = 125.04452 - 1934.136261 * T;

  // Nutations (arcsec)
  double dPSI = -17.20 * sind(OMEGA) - 1.32 * sind(2.0 * L) - 0.23 * sind(2.0 * dL) + 0.21 * sind(2.0 * OMEGA);
  double dEPSILON = 9.20 * cosd(OMEGA) + 0.57 * cosd(2.0 * L) + 0.10 * cosd(2.0 * dL) - 0.09 * cosd(2.0 * OMEGA);

  // arcsec -> deg
  dPSI /= 3600.0;
  dEPSILON /= 3600.0;

  const double gast = wrap360(THETAm + dPSI * cosd(EPSILONm + dEPSILON));
  return gast;
}

// Julian date of 2015-01-01 00:00:00 UTC
static constexpr double kJD2015_0101_0000 = 2457023.5;

double daysSince2015(double jd) {
  return jd - kJD2015_0101_0000;
}

}  // namespace geomag_demo
