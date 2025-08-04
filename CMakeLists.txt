cmake_minimum_required(VERSION 3.17)
project(cpp_FEM_multiphysics)

set(CMAKE_CXX_STANDARD 17)
include_directories(include)

# --- 依赖库管理 ---
# CMake 将首先检查 'share' 目录下是否存在所需的库。
# 如果不存在，它将自动从远程仓库下载，但不会创建 subbuild 子项目。

# --- Eigen ---
set(EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/share/eigen)
if(NOT EXISTS ${EIGEN3_INCLUDE_DIR})
        message(STATUS "本地未找到 Eigen，正在从远程下载...")
        execute_process(
                COMMAND git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git ${EIGEN3_INCLUDE_DIR}
                RESULT_VARIABLE GIT_CLONE_RESULT
        )
        if(NOT GIT_CLONE_RESULT EQUAL 0)
                message(FATAL_ERROR "克隆 Eigen 仓库失败")
        endif()
        message(STATUS "Eigen 下载完成")
endif()

# --- Boost ---
# 注意: Boost 仓库非常大，这里只克隆最新版本以减小体积
set(BOOST_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/share/boost)
if(NOT EXISTS ${BOOST_INCLUDE_DIR})
        message(STATUS "本地未找到 Boost，正在从远程下载...")
        execute_process(
                COMMAND git clone --depth 1 --branch boost-1.85.0 https://github.com/boostorg/boost.git ${BOOST_INCLUDE_DIR}
                RESULT_VARIABLE GIT_CLONE_RESULT
        )
        if(NOT GIT_CLONE_RESULT EQUAL 0)
                message(FATAL_ERROR "克隆 Boost 仓库失败")
        endif()
        message(STATUS "Boost 下载完成")
endif()

# --- OpenMP ---
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ot")


# --- 源文件 ---
file(GLOB_RECURSE SOURCES
        ${PROJECT_SOURCE_DIR}/src/core/*.cpp
        ${PROJECT_SOURCE_DIR}/src/physics/*.cpp
        ${PROJECT_SOURCE_DIR}/src/io/*.cpp
        ${PROJECT_SOURCE_DIR}/src/solver/*.cpp
        ${PROJECT_SOURCE_DIR}/src/utils/*.cpp
        ${PROJECT_SOURCE_DIR}/src/post/*.cpp
)

# --- 主库 (femlib) ---
add_library(femlib
        ${SOURCES}
)
target_link_libraries(femlib PUBLIC
        OpenMP::OpenMP_CXX
)
target_include_directories(femlib PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
        ${BOOST_INCLUDE_DIR}
)


# --- 主程序 ---
add_executable(cpp_FEM_multiphysics src/main.cpp)
target_link_libraries(cpp_FEM_multiphysics femlib)


# --- Google Test (推荐继续使用 FetchContent) ---
# Google Test 的集成比较复杂，FetchContent 能完美处理其编译和链接目标，
# 且其体积小，对加载速度影响微乎其微。
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


# --- 测试 ---
enable_testing()
include(GoogleTest)

# 创建测试可执行文件
file(GLOB_RECURSE TEST_SOURCES tests/*.cpp)
add_executable(run_tests ${TEST_SOURCES})

# 链接测试库
target_link_libraries(run_tests femlib gtest_main)
target_include_directories(run_tests PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIR}
)

# 自动发现测试
gtest_discover_tests(run_tests)


# --- 编译器特定设置 ---
if(MSVC)
        # 禁用 C4819 警告
        target_compile_options(femlib PRIVATE /wd4819)
        target_compile_options(cpp_FEM_multiphysics PRIVATE /wd4819)
        target_compile_options(run_tests PRIVATE /wd4819)
endif()