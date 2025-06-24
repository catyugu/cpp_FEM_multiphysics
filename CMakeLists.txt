cmake_minimum_required(VERSION 3.17)
project(cpp_FEM_multiphysics)
enable_language(Fortran)
include_directories(include)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Od" )
# Add the path to your Eigen installation to CMAKE_PREFIX_PATH
# This helps CMake find the Eigen3 package.
# If you installed to the default location, you might not need this.
include(FetchContent)
FetchContent_Declare(eigen GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0)
FetchContent_MakeAvailable(eigen)
include_directories(${eigen_SOURCE_DIR})

file(GLOB_RECURSE SOURCES
        ${PROJECT_SOURCE_DIR}/src/core/*.cpp
        ${PROJECT_SOURCE_DIR}/src/physics/*.cpp
        ${PROJECT_SOURCE_DIR}/src/io/*.cpp
        )
add_library(femlib
        ${SOURCES}
)
target_include_directories(femlib PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)
target_include_directories(femlib PUBLIC
        include
)

# --- Main Executable ---
add_executable(cpp_FEM_multiphysics src/main.cpp)
target_link_libraries(cpp_FEM_multiphysics femlib)

# 使用 FetchContent 引入 Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# --- Tests ---
enable_testing()
include(GoogleTest)

# Create the test executable
file(GLOB_RECURSE TEST_SOURCES tests/*.cpp)
add_executable(run_tests ${TEST_SOURCES})

# Link the test executable against our solver library and Google Test
target_link_libraries(run_tests femlib gtest_main)

# Automatically discover and add tests to the CTest runner
gtest_discover_tests(run_tests)