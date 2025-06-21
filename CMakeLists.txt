cmake_minimum_required(VERSION 3.17)
project(cpp_FEM_multiphysic)

set(CMAKE_CXX_STANDARD 17)

# Add the path to your Eigen installation to CMAKE_PREFIX_PATH
# This helps CMake find the Eigen3 package.
# If you installed to the default location, you might not need this.
set(EIGEN3_INCLUDE_DIR "E:/env/cpp/eigen-3.4.0")
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(cpp_FEM_multiphysic ${SOURCES})
file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/src/*.h ${PROJECT_SOURCE_DIR}/src/*.hpp)
target_include_directories(cpp_FEM_multiphysic PRIVATE ${PROJECT_SOURCE_DIR}/src)
# Link Eigen to your executable
target_include_directories(cpp_FEM_multiphysic PRIVATE ${EIGEN3_INCLUDE_DIR})


# 使用 FetchContent 引入 Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
# 添加测试可执行文件
add_executable(google_test test/test.cpp)
target_include_directories(
        google_test
        PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${EIGEN3_INCLUDE_DIR}
)
# 链接 Google Test 库
target_link_libraries(google_test PRIVATE gtest gtest_main)

# 启用测试
enable_testing()
add_test(NAME google_test COMMAND google_test)

