cmake_minimum_required(VERSION 3.17)
project(cpp_FEM_multiphysics)
include_directories(include)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED) # Find OpenMP
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ot")
set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/share/eigen-3.4.0)
set(BOOST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/share/boost-1.88.0)

if(EXISTS ${EIGEN3_INCLUDE_DIR})
        message(STATUS "Found local Eigen library at: ${EIGEN3_INCLUDE_DIR}")
else()
        message(FATAL_ERROR "Eigen library not found at the specified path: ${EIGEN3_INCLUDE_DIR}. Please download and place it there.")
endif()

# Check if the directory exists
if(EXISTS ${BOOST_INCLUDE_DIR})
        message(STATUS "Found local Boost library at: ${BOOST_INCLUDE_DIR}")
else()
        message(FATAL_ERROR "Boost library not found at the specified path: ${BOOST_INCLUDE_DIR}. Please download and place it there.")
endif()


file(GLOB_RECURSE SOURCES
        ${PROJECT_SOURCE_DIR}/src/core/*.cpp
        ${PROJECT_SOURCE_DIR}/src/physics/*.cpp
        ${PROJECT_SOURCE_DIR}/src/io/*.cpp
        ${PROJECT_SOURCE_DIR}/src/solver/*.cpp
        ${PROJECT_SOURCE_DIR}/src/utils/*.cpp
        ${PROJECT_SOURCE_DIR}/src/post/*.cpp
        )
add_library(femlib
        ${SOURCES}
)
target_link_libraries(femlib PUBLIC
        OpenMP::OpenMP_CXX
)


target_include_directories(femlib PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)
target_include_directories(femlib PUBLIC
        include
        ${EIGEN3_INCLUDE_DIR}
        ${BOOST_INCLUDE_DIR}
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
# 添加测试可执行文件
file(GLOB_RECURSE TEST_SOURCES tests/*.cpp)
add_executable(google_test ${TEST_SOURCES})
target_include_directories(
        google_test
        PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
        ${BOOST_INCLUDE_DIR}
)
# 链接 Google Test 库
target_link_libraries(google_test PRIVATE gtest gtest_main)

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

if(MSVC)
        # Disables warning C4819: The file contains a character that cannot be
        # represented in the current code page. Save the file in Unicode format to prevent data loss.
        target_compile_options(femlib PRIVATE /wd4819)
        target_compile_options(cpp_FEM_multiphysics PRIVATE /wd4819)
        target_compile_options(google_test PRIVATE /wd4819)
        target_compile_options(run_tests PRIVATE /wd4819)
endif()
