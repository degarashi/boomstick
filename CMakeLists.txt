cmake_minimum_required(VERSION 2.8)
project(boomstick)
include("crosscompile/${BUILD_TYPE}.cmake" require)
include_directories(${PROJECT_SOURCE_DIR}/spinner)
link_directories(${CMAKE_INSTALL_PREFIX}/lib_${TOOL_PREFIX})

file(GLOB MODEL3D "./*3D.cpp")
file(GLOB MODEL2D "./*2D.cpp")
add_library(boomstick STATIC ${MODEL2D} ${MODEL3D} model.cpp)
install(TARGETS boomstick ARCHIVE DESTINATION lib_${TOOL_PREFIX})
file(GLOB HEADERS "./*.hpp")
install(FILES ${HEADERS} DESTINATION include/boomstick)

if(NOT DEFINED ARCHITECTURE)
	# テストルーチンについては幾つかの警告を抑制する
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-unused-variable -Wno-unused-function")
	aux_source_directory(tests TEST_SRC)
	add_executable(bs_test spinner/main.cpp ${TEST_SRC})
	set(LINK_LIBS
		boomstick
		spinner
		boost_system
		boost_regex
		gtest
		pthread)
	target_link_libraries(bs_test ${LINK_LIBS})
	install(TARGETS bs_test RUNTIME DESTINATION bin_${TOOL_PREFIX})
endif()
