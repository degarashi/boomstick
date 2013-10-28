cmake_minimum_required(VERSION 2.8)
project(boomstick)

set(SSE_LEVEL 2)
include(spinner/CMakePC.txt require)
include(spinner/CMakeCompiler.txt require)
include_directories(${PROJECT_SOURCE_DIR}/spinner)
link_directories(${CMAKE_INSTALL_PREFIX}/lib_${TOOL_PREFIX})

add_library(boomstick STATIC geom2D.cpp col2D.cpp rigid2D.cpp collision.cpp)
add_executable(bs_test main.cpp)
target_link_libraries(bs_test boomstick spinner boost_system)

install(TARGETS boomstick bs_test
	ARCHIVE DESTINATION lib_${TOOL_PREFIX}
	RUNTIME DESTINATION bin_${TOOL_PREFIX})
file(GLOB HEADERS "./*.hpp")
install(FILES ${HEADERS} DESTINATION include/boomstick)
