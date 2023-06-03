cmake_minimum_required(VERSION 3.3)
project(lbfgspp VERSION 0.2.0 LANGUAGES CXX)

# + ----------------- +
# | BUILDING SETTINGS |
# + ----------------- +

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# + ----------------- +
# | COMPILATION FLAGS |
# + ----------------- +

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wall WALL_SUPPORTED)
if(${WALL_SUPPORTED})
  add_compile_options(-Wall)
endif()

# + --------------- +
# | LBFGSpp LIBRARY |
# + --------------- +

add_library(lbfgspp INTERFACE)

# + -------- +
# | INCLUDES |
# + -------- +

target_include_directories(lbfgspp INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# + ----------------------- +
# | FIND EXTERNAL LIBRARIES |
# + ----------------------- +

find_package(Eigen3 3.0 REQUIRED)
target_link_libraries(lbfgspp INTERFACE Eigen3::Eigen)
message("-- Eigen3 version: " ${EIGEN3_VERSION_STRING})

# + ------------ +
# | INSTALLATION |
# + ------------ +
