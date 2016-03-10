cmake_minimum_required(VERSION 2.8)
project(boomstick)
include("crosscompile/${BUILD_TYPE}.cmake" require)
include_directories(${PROJECT_SOURCE_DIR}/spinner)
link_directories(${CMAKE_INSTALL_PREFIX}/lib_${TOOL_PREFIX})

find_package(Boost REQUIRED system regex)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

file(GLOB_RECURSE HEADERS "*.hpp")
file(GLOB MODEL3D "./*3D.cpp")
file(GLOB MODEL2D "./*2D.cpp")
aux_source_directory(ntree NTREE_SRC)
set(LIB_SRC bits.cpp)
add_library(boomstick STATIC ${MODEL2D} ${MODEL3D} ${NTREE_SRC} ${LIB_SRC} model.cpp)
install(TARGETS boomstick ARCHIVE DESTINATION lib_${TOOL_PREFIX})

set(COPYHEADERS_EXCLUDE_DIR "spinner")		# ここに除外したいヘッダディレクトリを指定
set(COPYHEADERS_MY_DIR "boomstick")
include(crosscompile/CopyHeaders.cmake require)
include_directories(
	${Boost_INCLUDE_DIRS}
	${GTEST_INCLUDE_DIRS}
)

if(NOT DEFINED ARCHITECTURE)
	# テストルーチンについては幾つかの警告を抑制する
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-unused-variable -Wno-unused-function")
	aux_source_directory(tests TEST_SRC)
	add_executable(bs_test spinner/main.cpp ${TEST_SRC})
	set(LINK_LIBS
		boomstick
		spinner
		${Boost_LIBRARIES}
		${GTEST_LIBRARIES}
		${CMAKE_THREAD_LIBS_INIT}
	)
	target_link_libraries(bs_test ${LINK_LIBS})
	install(TARGETS bs_test RUNTIME DESTINATION bin_${TOOL_PREFIX})
endif()
