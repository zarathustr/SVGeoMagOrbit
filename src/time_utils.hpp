#pragma once

#include <string>

namespace geomag_demo {

struct DateTime {
  int year = 0;
  int month = 0;   // 1-12
  int day = 0;     // 1-31
  int hour = 0;
  int minute = 0;
  double second = 0.0;
};

// Parse an STK timestamp, typically like:
//   "01 Mar 2019 00:00:00.000"
// Returns true on success.
bool parseStkTimestamp(const std::string& ts, DateTime& out);

// Convert calendar date/time (UTC) to Julian date.
double toJulianDateUTC(const DateTime& t);

// MATLAB-compatible GMST and GAST, in degrees.
double jdToGMST_deg(double jd);
double jdToGAST_deg(double jd);

// Utility: days since 2015-01-01 00:00:00 UTC (IGRF 2015 epoch in this demo).
double daysSince2015(double jd);

}  // namespace geomag_demo
