# C++ FEM Multiphysics Solver

`cpp-fem-multiphysics` is a lightweight, header-only finite element analysis (FEA) framework built in modern C++. It is designed to solve a variety of physics problems, with a focus on multiphysics coupling, such as electro-thermal analysis.

-----

## Key Features

* **Finite Element Core**: A robust and extensible core architecture featuring classes for `Mesh`, `Node`, `Element`, and `DOFManager`.
* **Multiphysics Support**:
    * **Heat Transfer**: Solves 1D, 2D, and 3D steady-state and transient heat conduction problems.
    * **Electromagnetics**: Simulates 1D and 2D current distribution (Voltage).
    * **Coupled Physics**: Handles strongly coupled electro-thermal problems where electrical conductivity is temperature-dependent and Joule heating acts as a thermal source.
* **Modern C++**: Built using C++17 features for clean, efficient, and maintainable code.
* **Advanced Solvers**: Utilizes the powerful **Eigen** library for sparse linear algebra and provides a `SolverFactory` to select the appropriate solution strategy (single-field vs. coupled).
* **Flexible Boundary Conditions**: Supports Dirichlet (fixed value), Neumann (flux), and Cauchy (mixed/convection) boundary conditions.
* **Mesh Handling**: Includes built-in generators for uniform 1D, 2D, and 3D meshes and an importer for COMSOL's `.mphtxt` file format.
* **Testing**: A comprehensive test suite using GoogleTest verifies the correctness of individual components and coupled simulations.
* **Exporting**: Exports results to `.vtk` files for visualization in standard tools like ParaView or VisIt.

-----

## Project Structure

The project is organized into a clean and logical directory structure:

```
cpp_FEM_multiphysics/
├── CMakeLists.txt          # Main build script
├── README.md               # You are here!
├── include/                # Header files (.hpp)
│   ├── core/               # Core architectural components (Mesh, Problem, Element)
│   ├── io/                 # Input/Output utilities (VTK Exporter, COMSOL Importer)
│   ├── physics/            # Physics-specific modules (Heat2D, Current1D, etc.)
│   ├── solver/             # Solution strategy classes (SingleFieldSolver, CoupledSolver)
│   └── utils/              # General utilities (SimpleLogger)
├── src/                    # Source files (.cpp)
│   ├── core/
│   ├── io/
│   ├── physics/
│   ├── solver/
│   └── main.cpp            # Main application entry point
├── tests/                  # Google Test source files for validation
└── docs/                   # Detailed markdown documentation for each namespace
```

-----

## Requirements

* C++17 compatible compiler (e.g., GCC, Clang, MSVC)
* CMake \>= 3.17
* Eigen \>= 3.3
* GoogleTest (will be fetched automatically by CMake)

-----

## Getting Started

### Build Instructions

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd cpp_FEM_multiphysics
    ```
2.  **Configure Eigen:**
    Open `CMakeLists.txt` and update the `EIGEN3_INCLUDE_DIR` variable to point to the root directory of your Eigen library installation.
    ```cmake
    # On line 9
    set(EIGEN3_INCLUDE_DIR "path/to/your/eigen-3.4.0")
    ```
3.  **Build the project:**
    ```bash
    mkdir build && cd build
    cmake ..
    make
    ```
    This will create two main executables: `cpp_FEM_multiphysics` and `run_tests`.

### Running the Simulation

The main executable runs a pre-configured simulation of a 2D heat conduction problem on a circular mesh imported from a COMSOL file.

To run it:

```bash
./cpp_FEM_multiphysics
```

This will generate a log file `femsolver.log` and a result file `comsol_circle_steadystate_results.vtk` that you can open in ParaView.

### Running Tests

The project includes a full suite of tests to validate the functionality of the solvers.

To run the tests:

```bash
./run_tests
```

All tests should pass, confirming that the FEM implementation is correct against analytical solutions and benchmarks.

-----

## Detailed Documentation

For a deeper dive into the architecture and API, please refer to the detailed documentation for each namespace:

* **[Core Namespace](docs/core.md)**: Describes the fundamental building blocks of the FEM framework.
* **[PhysicsField Namespace](docs/physics.md)**: Details the specific implementations for different physical simulations.
* **[Solver Namespace](docs/solver.md)**: Explains the available solution strategies.
* **[IO Namespace](docs/io.md)**: Covers the utilities for importing and exporting data.