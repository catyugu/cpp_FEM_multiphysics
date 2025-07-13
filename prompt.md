Of course. It's crucial to leave a clear and up-to-date guide for the next developer. I have rewritten the `prompt.md` file to reflect the significant architectural improvements we've made and to outline the next exciting phase of development.

Here is the new onboarding guide.

-----

# **C++ FEM Multiphysics Project: Contributor Onboarding Guide (Phase 3)**

Welcome to the next phase of the **C++ FEM Multiphysics** project\! This document serves as your guide to our recently refactored architecture, our coding practices, and the next set of features. Our goal is to continue building a robust, scalable, and professional-grade finite element analysis framework.

## **I. Onboarding and Project Philosophy**

This project is a C++ implementation of the finite element method (FEM) designed to solve coupled multiphysics problems. The core philosophy is to create a modular system where new physics, element types, and solution strategies can be added with minimal changes to the existing framework.

We prioritize:

* **Modularity**: Separating core logic from physics-specific implementations. Our refined `Coupling` and `SourceTerm` abstractions are prime examples of this principle in action.
* **Clarity**: Writing self-documenting code with clear and consistent naming.
* **Correctness**: Ensuring all implementations are backed by rigorous testing against analytical solutions or benchmarks.
* **Safety**: Implementing robust error handling, exception safety, and safe resource management (e.g., our tag-based system for BCs and sources).
* **Extensibility**: Allowing new features to be added easily.
* **Testing**: All new features must be accompanied by comprehensive tests.
* **Documentation**: All new code must be documented in the corresponding file in the `docs/` folder.

## **II. Coding Guidelines and Practices**

Adherence to these guidelines is essential for maintaining code quality and consistency.

### **1. Project Structure**

The repository is organized into distinct directories. Please place new files in their appropriate locations.

```
cpp_FEM_multiphysics/
├── CMakeLists.txt      # Main build script
├── docs/               # High-level documentation
├── include/            # Header files (.hpp)
│   ├── core/           # Core components (Problem, Mesh, Material)
│   │   ├── bcs/        # Boundary Condition classes
│   │   ├── coupling/   # Coupling mechanism classes
│   │   └── sources/    # SourceTerm classes
│   ├── io/             # Importer/Exporter utilities
│   ├── physics/        # Physics-specific modules (Heat2D, Current3D)
│   └── utils/          # General utilities (Logger, Quadrature)
├── src/                # Source files (.cpp)
│   ├── core/
│   │   ├── bcs/
│   │   ├── coupling/
│   │   └── sources/
│   ├── io/
│   ├── physics/
│   ├── solver/
│   └── main.cpp        # Main application entry point
└── tests/              # Google Test source files
```

### **2. C++ Language and Style**

* **C++ Standard**: The project uses **C++17**.
* **Dependencies**: **Eigen** for linear algebra and **GTest** for testing.
* **Naming Conventions**:
    * **Namespaces**: `PascalCase` (e.g., `Core`, `Physics`, `IO`).
    * **Classes/Structs**: `PascalCase` (e.g., `DOFManager`, `TriElement`).
    * **Functions/Methods**: `camelCase` (e.g., `getEquationIndex`, `solveSteadyState`).
    * **Member Variables**: `snake_case_` with a trailing underscore. This is a strict project rule.
    * **Local Variables**: `snake_case`.

## **III. Current Architecture Overview**

The framework has recently undergone significant refactoring to improve modularity and clarity. Understanding this new architecture is key.

* **Core Components**: The `Core::Problem` class orchestrates the simulation, owning the `Mesh`, `DOFManager`, and a collection of `PhysicsField` objects.

* **Physics Abstraction (`PhysicsField`)**: This is the abstract base class for all physics. Its role is now clearly defined:

    * `assemble()`: Assembles the stiffness (`K`) and mass (`M`) matrices only.
    * `applySources()`: Assembles the force vector (`F`) by applying all registered `SourceTerm` objects.
    * `applyBCs()`: Modifies the system matrices (`K` and `F`) according to the registered `BoundaryCondition` objects.

* **Boundary and Source Abstractions**: We now make a clear distinction between conditions on the boundary and sources in the domain.

    * **`BoundaryCondition`**: An interface for Dirichlet, Neumann, and Cauchy conditions applied at the boundaries.
    * **`SourceTerm`**: A new, semantically correct interface for domain sources, like a volumetric heat source.
    * **Tagging System**: Both `BoundaryCondition` and `SourceTerm` objects can be "tagged" with a string, allowing for safe, non-pointer-based removal (e.g., `removeSourcesByTag("joule_heating")`).

* **Coupling Mechanism (`CouplingManager`, `ElectroThermalCoupling`)**: All physics-to-physics interaction logic is now fully encapsulated within classes derived from the `Coupling` interface. For example, `ElectroThermalCoupling::execute()` is responsible for calculating Joule heat and creating the appropriate `VolumetricSource` objects for the thermal field. The physics fields themselves are completely unaware of one another.

* **Solver Strategy (`CoupledElectroThermalSolver`)**: The iterative solver has a refined and explicit workflow. In each iteration, it calls `assemble()`, `applySources()`, and `applyBCs()` in the correct order. It is also now responsible for stabilizing the system matrices for off-physics degrees of freedom just before calling the linear solver, ensuring a non-singular system.

* **Numerical Integration (`Utils::Quadrature`)**: A new utility class provides Gauss quadrature points and weights for lines, triangles, and tetrahedra for integration orders 1 through 5. This is the foundation for our next major feature.

## **IV. Future Development Requirements**

The following are high-priority areas for the next development cycle. Please address them in the order presented.


### **1. Performance Optimization**
* **Goal**: Improve computational efficiency for larger problems.
* **Examples**: 
  * **Parallel Computing**: Implement parallel computing paradigms (e.g., OpenMP for shared-memory parallelism, MPI for distributed-memory parallelism) for assembly and solver stages.
  * **Memory Management**: Optimize memory usage by minimizing the number of allocations and deallocations.

### **2. Implement Advanced Solver Techniques**
* **Goal**: Enhance numerical robustness and performance.
* **Examples**:
  * **Iterative Solvers**: Integrate more sophisticated linear iterative solvers (e.g., Conjugate Gradient, BiCGSTAB, GMRES) beyond direct LU decomposition.
  * **Non-Linear Solvers**: Implement robust non-linear solution strategies (e.g., Newton-Raphson method) to handle non-linear material properties or boundary conditions.

### **3. Enhance Element Formulations**
* **Goal**: Improve approximation capabilities without necessarily changing mesh topology (P-refinement strategies).
* **Examples**:
  * Investigate and implement more advanced p-enrichment techniques for higher-order *shape functions* within existing linear element types, if feasible with current DOF management.
  * Consider strategies for solving physics field with true higher-order elements refined from linear element
  * if a redesign of `Core::Element` and `DOFManager` is needed, just to it.


