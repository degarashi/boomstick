cmake_minimum_required(VERSION 2.8)
project(boomstick)
set(CMAKE_CXX_COMPILER "i386-mingw32-g++")
set(CMAKE_AR "i386-mingw32-ar")
set(CMAKE_RANLIB "i386-mingw32-ranlib")
include_directories("/home/slice/local/include/"
					${PROJECT_SOURCE_DIR}
					${PROJECT_SOURCE_DIR}/spinner)
link_directories(/tmp/spinner_build/)
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -std=c++11 -ggdb3 -O0 -msse2 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -msse2")
set(CMAKE_BUILD_TYPE "Release")
add_library(boomstick STATIC geom2D.cpp col2D.cpp rigid2D.cpp collision.cpp)
