#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace geomag_demo {

struct StkRow {
  std::string timestamp;              // first 24 chars (STK "Time (UTCG)" column)
  std::vector<double> values;         // numeric columns
};

// Load an STK "Report" style text file.
// The parser searches for a header line containing "Time (UTCG)", skips one line,
// then expects each data line to start with a fixed-width timestamp field.
//
// Parameters
// - file: path to the STK report
// - max_rows: maximum number of data lines to load (0 => load all)
// Returns true on success.
bool loadStkReport(const std::filesystem::path& file,
                   std::vector<StkRow>& out,
                   std::size_t max_rows = 0);

}  // namespace geomag_demo
