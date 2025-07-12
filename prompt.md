# C++ FEM Multiphysics Project: Contributor Onboarding Guide

Welcome to the **C++ FEM Multiphysics** project\! This document serves as your guide to understanding the existing architecture, adhering to our coding practices, and tackling the next set of complex features. Our goal is to build a robust, scalable, and maintainable finite element analysis framework. Please read this guide carefully before writing any code.

## Onboarding and Project Philosophy

This project is a C++ implementation of the finite element method (FEM) designed to solve coupled multiphysics problems. The core philosophy is to create a modular system where new physics, element types, and solution strategies can be added with minimal changes to the existing framework.

We prioritize:

* **Modularity**: Separating core logic from physics-specific implementations.
* **Clarity**: Writing self-documenting code with clear and consistent naming.
* **Correctness**: Ensuring all implementations are backed by rigorous testing against analytical solutions or benchmarks.
* **Safety** : Implementing robust error handling and exception safety.
* **Performance**: Optimizing for speed and memory efficiency.
* **Extensibility**: Allowing new physics, element types, and solution strategies to be added easily.
* **Testing**: All new features must be accompanied by comprehensive tests.
* **Documentation**: All new code must be documented in the corresponding file in the folder docs/ .
-----

## I. Coding Guidelines and Practices

Adherence to these guidelines is essential for maintaining code quality and consistency.

### 1\. Project Structure

The repository is organized into distinct directories. Please place new files in their appropriate locations.

```
cpp_FEM_multiphysics/
├── CMakeLists.txt      # Main build script
├── docs/               # High-level documentation (core.md, physics.md)
├── include/            # Header files (.hpp)
│   ├── core/           # Core architectural components (Problem, Mesh, Solver)
│   ├── io/             # Importer/Exporter utilities
│   ├── physics/        # Physics-specific modules (Heat2D, Current1D)
│   └── utils/          # General utilities (SimpleLogger)
├── src/                # Source files (.cpp)
│   ├── core/
│   ├── io/
│   ├── physics/
│   └── main.cpp        # Main application entry point
└── tests/              # Google Test source files
```

### 2\. C++ Language and Style

* **C++ Standard**: The project uses **C++17**. Use C++17 features where they improve clarity and safety.
* **Dependencies**: The primary dependencies are **Eigen** for linear algebra and **GTest** for testing.
* **Naming Conventions**:
    * **Namespaces**: `PascalCase` (e.g., `Core`, `Physics`, `IO`). All code should reside within a namespace.
    * **Classes/Structs**: `PascalCase` (e.g., `DOFManager`, `TriElement`).
    * **Functions/Methods**: `camelCase` (e.g., `getEquationIndex`, `solveSteadyState`).
    * **Member Variables**: `snake_case_` with a trailing underscore (e.g., `mesh_`, `num_equations_`). This is a strict project rule.
    * **Local Variables**: `snake_case` (e.g., `k_triplets`, `num_steps`).

### 3\. Core Architectural Principles

* **Problem Class**: The `Core::Problem` class is the main orchestrator. It owns the mesh, DOF manager, and physics fields. It does **not** contain direct solving logic.
* **Solver Strategy**: The actual solution logic is delegated to `Solver` objects (e.g., `SingleFieldSolver`, `CoupledElectroThermalSolver`). The correct solver is chosen by the `SolverFactory`. When adding new coupled physics, create a new `Solver` strategy.
* **PhysicsField Abstraction**: All physics (e.g., thermal, electrical) must inherit from the `Physics::PhysicsField` base class. This class defines the interface for assembling matrices (`K`, `M`) and the RHS vector/matrix (`F`).
* **Memory Management**:
    * The `Problem` class uses `std::unique_ptr` to manage the lifetime of the `Mesh`, `DOFManager`, and `PhysicsField` objects.
    * The `Mesh` class uses raw pointers for nodes and elements, as it acts as a non-owning container. The `Problem`'s destructor is responsible for their cleanup. Follow this pattern carefully to avoid memory leaks.
* **Matrix/Vector Operations**: All linear algebra must use the **Eigen** library. The RHS (`F_`) and solution (`U_`) are `Eigen::MatrixXd` to support multiple load cases.

### 4\. Testing

* **Framework**: All tests must be written using **Google Test**.
* **Location**: Place new test files in the `/tests` directory.
* **Practice**: For any new feature or bug fix, a corresponding test must be created. Validate new physics implementations against known analytical solutions (as seen in `test_heat_3D.cpp` and `test_current_2D.cpp`) or expected physical behavior (as in `test_coupled_2D.cpp`).

### 5\. Documentation and Logging

* **High-Level Docs**: Major new components should be documented in the `/docs` directory in Markdown format.
* **Code Comments**: Use Doxygen-style comments for headers to explain class and method purposes.
* **Logging**: Use the provided `SimpleLogger` for all console output. Avoid using `std::cout`. Log key steps in the setup, assembly, and solve processes.

-----

## II. Future Development Requirements

The following are high-priority areas for future development.

### 1\. Advanced Mesh IO

* **Problem**: The `Importer` and `Exporter` only supports a specific COMSOL text format and scalar output.
* **Requirement**: Enhance the framework to support more complex geometries.
* **Tasks**:
    1.  **Importer Upgrade**: Edit the class `Importer` that can read `.msh` files from the popular open-source mesher Gmsh. This will involve parsing its file format for nodes, elements (lines, triangles, tets), and physical groups (for assigning boundaries).
    2.  **Exporter Upgrade**: Edit the class `Exporter` that can write `.vtk` files for visualization. This will involve writing the node coordinates and element connectivity.

### 2\. New Physics: 2D Magnetic

* **Problem**: The framework currently handles thermal and electromagnetic physics.
* **Requirement**: Add a new physics module for **2D Magnetic field**.
* **Tasks**:
    1.  **New PhysicsField**: Create `Magnetic2D`class inheriting from `Physics::PhysicsField`. The primary variable will be "Displacement," which is a vector quantity (u, v, w).
    2.  **DOFManager Update**: The `DOFManager` currently assumes one variable per node (e.g., "Temperature"). It must be updated to handle vector variables like displacement, registering "Displacement\_X", "Displacement\_Y", etc., for each node.
    3.  **New Material Properties**: The `Material` class needs to support properties needed for the .
    4.  **New Element Formulation**: The `assemble` method in the new solid mechanics classes will compute the element stiffness matrix `k_e = integral(B^T * D * B) dV`, where `B` is the strain-displacement matrix and `D` is the constitutive matrix (stress-strain relationship).

### 3\. Implement Adaptive Mesh Refinement (AMR)

* **Problem**: Implement a posteriori error estimation and an adaptive mesh refinement loop to automatically improve solution accuracy in regions with high error gradients.

* **Tasks**:
    1.  **Error Estimator**: Implement a **Zienkiewicz-Zhu error estimator**. This involves:
        * First, solving the primary problem for the solution (e.g., temperature).
        * Then, projecting the discontinuous gradients (stresses/fluxes) from the integration points to the nodes to create a smoothed, continuous gradient field.
        * Calculating an error norm for each element based on the difference between the smoothed nodal gradients and the original discontinuous element gradients.
    2.  **Refinement Strategy**:
        * Based on the error estimates, mark a certain percentage of elements (e.g., the top 25% with the highest error) for refinement.
    3.  **Mesh Management**:
        * The `Core::Problem` class must be modified to handle a mesh that changes during the solution. The `DOFManager` will need to be rebuilt after each refinement step.
        * You will need to implement a function to **transfer the solution** from the old mesh to the newly refined mesh to provide a good initial guess for the next solve step.
    4.  **Validation Test**: Create a test case for a classic AMR benchmark: a 2D L-shaped domain heat transfer problem. The analytical solution has a singularity in the re-entrant corner. Your AMR implementation should demonstrate that elements are selectively refined near this corner, and the solution should converge to the analytical value faster than with uniform mesh refinement.
