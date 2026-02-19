cmake_minimum_required(VERSION 3.16)

project(geomag_sv_demo LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# User option: path to directory containing the STK "Satellite2*.txt" files.
# You can set it at configure time:
#   cmake -S . -B build -DSAT_DATA_DIR=/path/to/stk_exports
set(SAT_DATA_DIR "" CACHE PATH "Directory containing STK Satellite2*.txt files (optional)")

# Eigen dependency
option(GEOMAG_DEMO_FETCH_EIGEN "Fetch Eigen automatically if not found" ON)
find_package(Eigen3 3.3 QUIET NO_MODULE)

if(NOT Eigen3_FOUND AND GEOMAG_DEMO_FETCH_EIGEN)
  include(FetchContent)
  message(STATUS "Eigen3 not found; fetching Eigen 3.4.0 via FetchContent")
  FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
  )
  FetchContent_MakeAvailable(eigen)
endif()

# Executable
add_executable(geomag_demo
  src/main.cpp
  src/stk_loader.cpp
  src/time_utils.cpp
  src/eci_ecef.cpp
  src/igrf.cpp
  src/orbit_propagator.cpp
)

target_include_directories(geomag_demo PRIVATE src)

if(TARGET Eigen3::Eigen)
  target_link_libraries(geomag_demo PRIVATE Eigen3::Eigen)
elseif(TARGET eigen)
  # Some Eigen CMake configs export target 'eigen' (rare). Keep as fallback.
  target_link_libraries(geomag_demo PRIVATE eigen)
else()
  message(FATAL_ERROR "Eigen not found. Install Eigen3 (e.g., libeigen3-dev) or set GEOMAG_DEMO_FETCH_EIGEN=ON")
endif()

# Provide a compile-time default data dir (overridable at runtime with --data_dir)
if(SAT_DATA_DIR STREQUAL "")
  target_compile_definitions(geomag_demo PRIVATE GEOMAG_DEMO_DATA_DIR="")
else()
  file(TO_CMAKE_PATH "${SAT_DATA_DIR}" SAT_DATA_DIR_CMAKE)
  target_compile_definitions(geomag_demo PRIVATE GEOMAG_DEMO_DATA_DIR="${SAT_DATA_DIR_CMAKE}")
endif()

# Copy coefficient files to build dir for convenience
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/igrfSg2015.txt ${CMAKE_CURRENT_BINARY_DIR}/igrfSg2015.txt COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/igrfSh2015.txt ${CMAKE_CURRENT_BINARY_DIR}/igrfSh2015.txt COPYONLY)

