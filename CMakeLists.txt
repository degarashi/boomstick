cmake_minimum_required(VERSION 2.8)
project(boomstick)

include_directories(${PROJECT_SOURCE_DIR}/spinner)
link_directories(/tmp/spinner_build/)
set(CMAKE_BUILD_TYPE "Debug")
include(spinner/CMake0.txt require)

add_library(boomstick STATIC geom2D.cpp col2D.cpp rigid2D.cpp collision.cpp)
add_executable(bs_test main.cpp)
target_link_libraries(bs_test boomstick spinner boost_system)
