cmake_minimum_required(VERSION 3.23)
project(arap)
# FetchContent added in CMake 3.11, downloads during the configure step
include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if (${CMAKE_VERSION} VERSION_LESS 3.14)
    include(cmake/add_FetchContent_MakeAvailable.cmake)
endif ()

FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
)

set(EIGEN_BUILD_DOC OFF)
# note: To disable eigen tests,
# you should put this code in a add_subdirectory to avoid to change
# BUILD_TESTING for your own project too since variables are directory
# scoped
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
FetchContent_MakeAvailable(eigen)

#set(eigen_include_dir ${CMAKE_BINARY_DIR}/eigen-src)
set(eigen_include_dir ${eigen_SOURCE_DIR})
set(CMAKE_CXX_STANDARD 11)
include_directories(include ${eigen_SOURCE_DIR})
add_executable(arap main.cpp src/Optimizer.cpp src/mesh.cpp include/common.h include/ARAP.h include/parameters.h src/quartic.cpp)
