#include "stk_loader.hpp"

#include <fstream>
#include <sstream>

namespace geomag_demo {

static bool containsTimeHeader(const std::string& line) {
  return line.find("Time (UTCG)") != std::string::npos;
}

static std::vector<double> parseDoubles(const std::string& s) {
  std::vector<double> v;
  std::istringstream iss(s);
  double x;
  while (iss >> x) {
    v.push_back(x);
  }
  return v;
}

bool loadStkReport(const std::filesystem::path& file,
                   std::vector<StkRow>& out,
                   std::size_t max_rows) {
  out.clear();

  std::ifstream ifs(file);
  if (!ifs.is_open()) {
    return false;
  }

  std::string line;
  bool found_header = false;

  // Find the header line containing "Time (UTCG)"
  while (std::getline(ifs, line)) {
    if (containsTimeHeader(line)) {
      found_header = true;
      break;
    }
  }

  if (!found_header) {
    return false;
  }

  // Skip the next line (usually a separator)
  std::getline(ifs, line);

  std::size_t count = 0;
  while (std::getline(ifs, line)) {
    if (max_rows != 0 && count >= max_rows) {
      break;
    }

    // Require at least timestamp field
    if (line.size() < 24) {
      continue;
    }

    StkRow row;
    row.timestamp = line.substr(0, 24);
    const std::string rest = line.substr(24);
    row.values = parseDoubles(rest);

    if (!row.values.empty()) {
      out.push_back(std::move(row));
      ++count;
    }
  }

  return !out.empty();
}

}  // namespace geomag_demo
