# C++ FEM Multiphysics Solver

`cpp-fem-multiphysics` is a lightweight, header-only finite element analysis (FEA) framework built in modern C++. It is designed to solve a variety of physics problems, with a strong focus on advanced features like p-refinement and multiphysics coupling, such as electro-thermal analysis.

-----

## Key Features

* **Finite Element Core**: A robust and extensible core architecture featuring classes for `Mesh`, `Node`, `Element`, `DOFManager`, `ElementGeometry`, and `FEValues`.
* **P-Refinement (Higher-Order Elements)**: Supports higher-order mathematical approximations (currently up to order 2 for triangles/tetrahedra, and order 5 for lines) on geometrically linear meshes, improving accuracy without remeshing.
* **Advanced DOF Management**: `DOFManager` intelligently handles degrees of freedom for both vertex and higher-order (edge) nodes for p-refinement.
* **Multiphysics Support**:
    * **Heat Transfer**: Solves 1D, 2D, and 3D steady-state and transient heat conduction problems.
    * **Electromagnetics**: Simulates 1D, 2D, and 3D current distribution (Voltage).
    * **Magnetostatics**: Solves 1D, 2D, and 3D magnetostatic problems using the magnetic potential.
    * **Coupled Physics**: Handles strongly coupled electro-thermal problems where electrical conductivity is temperature-dependent and Joule heating acts as a thermal source.
* **Modern C++**: Built using C++17 features for clean, efficient, and maintainable code.
* **Advanced Solvers**: Utilizes the powerful **Eigen** library for sparse linear algebra (supporting `SparseLU` and `BiCGSTAB` solvers) and provides a `SolverFactory` to select the appropriate solution strategy (single-field vs. coupled).
* **Flexible Boundary Conditions**: Supports Dirichlet (fixed value), Neumann (flux), and Cauchy (mixed/convection) boundary conditions, with robust handling for higher-order elements.
* **Material Properties**: Supports both constant and **temperature-dependent** material properties.
* **Mesh Handling**: Includes built-in generators for uniform 1D, 2D, and 3D meshes and importers for COMSOL's `.mphtxt` and **Gmsh's `.msh`** file formats.
* **Testing**: A comprehensive test suite using GoogleTest verifies the correctness of individual components and coupled simulations.
* **Exporting**: Exports results to `.vtk` files for visualization in standard tools like ParaView or VisIt, supporting both nodal (point) data and different element types (lines, triangles, and tetrahedra).

-----

## Project Structure

The project is organized into a clean and logical directory structure:

```

cpp\_FEM\_multiphysics/
├── CMakeLists.txt          \# Main build script
├── README.md               \# You are here\!
├── include/                \# Header files (.hpp)
│   ├── core/               \# Core architectural components (Mesh, Problem, Element)
│   │   ├── bcs/            \# Boundary Condition classes
│   │   ├── coupling/       \# Coupling mechanism classes
│   │   └── sources/        \# SourceTerm classes
│   ├── io/                 \# Input/Output utilities (VTK Exporter, COMSOL Importer)
│   ├── physics/            \# Physics-specific modules (Heat2D, Current1D, etc.)
│   ├── solver/             \# Solution strategy classes (SingleFieldSolver, CoupledSolver)
│   └── utils/              \# General utilities (SimpleLogger, Quadrature, ShapeFunctions, Exceptions)
├── src/                    \# Source files (.cpp)
│   ├── core/
│   │   ├── bcs/
│   │   ├── coupling/
│   │   └── sources/
│   ├── io/
│   ├── physics/
│   ├── solver/
│   └── main.cpp            \# Main application entry point
├── tests/                  \# Google Test source files for validation
└── docs/                   \# Detailed markdown documentation for each namespace

````

-----

## Requirements

* C++17 compatible compiler (e.g., GCC, Clang, MSVC)
* CMake >= 3.17
* Eigen >= 3.3
* GoogleTest (will be fetched automatically by CMake)

-----

## Getting Started

### Build Instructions

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd cpp_FEM_multiphysics
    ```

2.  **构建项目:**
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
````

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
* **[Core Mesh Components](docs/core/mesh.md)**: Details the mesh structure, elements, and the "Smart Element" refactoring with `ElementGeometry` and `FEValues`.
* **[Core Coupling Components](docs/core/coupling.md)**: Explains the generalized coupling mechanism.
* **[PhysicsField Namespace](docs/physics.md)**: Details the specific implementations for different physical simulations and the p-refinement strategy.
* **[Solver Namespace](docs/solver.md)**: Explains the available solution strategies and their enhancements.
* **[IO Namespace](docs/io.md)**: Covers the utilities for importing and exporting data, including updated VTK export capabilities.
* **[Utils Namespace](docs/utils.md)**: Provides general-purpose utilities like `ShapeFunctions`, `Quadrature`, and `SimpleLogger`.
* **[Exceptions Namespace](docs/exceptions.md)**: Describes the custom exception types for robust error handling.

-----

## Future Development

We are continually expanding the capabilities of this framework. Here are the key areas for future development:

### 1\. Implement Electromagnetic Field (Maxwell's Equations)

* **Goal**: Develop a full electromagnetic field simulation capability.
* **Focus**: The initial step is implementing a 3D Magnetostatics formulation (`Magnetostatics3D`) which will require upgrading the core components, especially the `DOFManager` and `PhysicsField` base class, to handle **vector-valued variables** (e.g., the magnetic vector potential with 3 components).

### 2\. Support for Advanced Research Problem Types (Solvers)

* **Goal**: Broaden the range of problems the solver can tackle for research applications.
* **Focus**: This includes implementing solvers for:
    * **Frequency Domain Steady-State**: For time-harmonic (AC) problems, requiring complex-number linear algebra.
    * Further enhancing the robustness and efficiency of existing **Coupled Transient Solvers**.

### 3\. Generalize Variable Material Properties for Coupled Fields

* **Goal**: Make material properties dynamically dependent on any relevant coupled field solution.
* **Current Status**: Temperature-dependent properties (e.g., electrical conductivity depending on temperature) are already implemented and used in electro-thermal coupling.
* **Future Focus**: Develop a more flexible mechanism within the `Material` class to evaluate properties based on a general set of input field values (e.g., allowing properties to depend on pressure, concentration, or other coupled variables).

-----