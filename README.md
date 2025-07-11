# C++ FEM Multiphysics
## Description
This project is a C++ implementation of a finite element method (FEM) for solving a multiphysics problem.
## Requirements
* C++17
* CMake >= 3.17
* Eigen >= 3.3
* GTest >= 1.10
## Build Guide
```bash
mkdir cmake-build-debug && cd cmake-build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
./run_tests ## Run google tests
./cpp_FEM_multiphysics ## Run main program
```
## Current Functions
* Finite Element Method Core.
* Solving 1D Heat&ElectroMagnetic Field and their coupling.
* Solving 2D Heat&ElectroMagnetic Field and their coupling.
* Export result in vtk file.
## Core File Structure
```bash
cpp_FEM_multiphysics
├── CMakeLists.txt      # Main build script
├── include/            # Header files (.hpp)
│   ├── core/           # Core architectural components
│   ├── io/             # Input/Output utilities (e.g., Exporter)
│   ├── physics/        # Physics-specific modules
│   └── utils/          # General utilities (e.g., Logger)
├── src/                # Source files (.cpp)
│   ├── core/
│   ├── io/
│   ├── physics/
│   └── main.cpp        # Main application entry point
└── tests/              # Google Test source files
```