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
```tree
cpp\_FEM\_multiphysics/
├── CMakeLists.txt      \# Main build script
├── docs/               \# High-level documentation
├── include/            \# Header files (.hpp)
│   ├── core/           \# Core components (Problem, Mesh, Material)
│   │   ├── bcs/        \# Boundary Condition classes
│   │   ├── coupling/   \# Coupling mechanism classes
│   │   └── sources/    \# SourceTerm classes
│   ├── io/             \# Importer/Exporter utilities
│   ├── physics/        \# Physics-specific modules (Heat2D, Current3D)
│   └── utils/          \# General utilities (Logger, Quadrature)
├── src/                \# Source files (.cpp)
│   ├── core/
│   │   ├── bcs/
│   │   ├── coupling/
│   │   └── sources/
│   ├── io/
│   ├── physics/
│   ├── solver/
│   └── main.cpp        \# Main application entry point
└── tests/              \# Google Test source files
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

The framework has evolved significantly. Understanding this new architecture, which enables higher-order approximations on simple linear meshes, is key.

### **1. P-Refinement Strategy**

The most significant recent upgrade is the implementation of a **p-refinement** (or p-enrichment) strategy.

* **Fixed Geometric Mesh**: The `Mesh` is always composed of simple, geometrically linear elements (e.g., a triangle always has 3 vertex nodes, a tetrahedron has 4). This decouples the geometric representation from the mathematical approximation.
* **Mathematical Order**: The accuracy of the simulation is controlled by the `element_order` property within each `PhysicsField`. Setting this to `2` or higher instructs the physics assembly to use more complex shape functions over the simple geometric element, effectively creating a higher-order approximation.
* **"Virtual" Nodes**: The `DOFManager` now intelligently creates and manages degrees of freedom for "virtual" nodes that don't exist in the mesh file but are mathematically required for higher-order elements. It currently supports nodes on the midpoints of element edges for 2nd-order approximations.

### **2. Advanced Degree of Freedom (DOF) Management**

To support p-refinement, the `DOFManager` has been completely redesigned.

* **Multiple DOF Maps**: It now maintains separate maps for different types of DOFs: one for vertices (`vertex_dof_map_`) and one for edges (`edge_dof_map_`).
* **Unique Edge Identification**: To uniquely identify a higher-order DOF on an edge, it uses a key composed of the **sorted** vertex IDs of that edge. This guarantees that `(Node1, Node2)` and `(Node2, Node1)` refer to the same edge DOF.
* **Order-Aware Build Process**: The `build()` method now requires the element orders of all physics fields to determine exactly which higher-order DOFs need to be created.

### **3. Robust Boundary Condition (BC) System**

Applying constraints to higher-order elements required a more sophisticated BC system.

* **Dual-Constructor `DirichletBC`**: The `DirichletBC` class now has two constructors: one that uses a `node_id` (for vertices) and a second that directly accepts an `equation_index`. This second constructor is essential for constraining the "virtual" higher-order DOFs, which do not have a node ID.
* **Consolidated Application**: The `PhysicsField::applyBCs()` method is now intelligent. It first consolidates all BCs into a map where each DOF index appears only once. This prevents common errors from redundantly applying constraints to the same node (e.g., a corner node shared by multiple boundary edges), which previously caused solver failures.

### **4. Solver and Performance**

* **Mixed-Order Coupled Solver**: The `CoupledElectroThermalSolver` has been upgraded to correctly stabilize and solve systems where different physics fields have different element orders.

---

## **III. Reconstruction Requirements**
The project's current structure provides an **excellent foundation** for a mature FEM core, demonstrating a strong grasp of separation of concerns and recent advancements like p-refinement. To evolve into a truly **mature and scalable FEM core**, further architectural shifts are necessary, primarily by making the `Element` class the central, intelligent entity in the assembly process.

The core philosophy for a mature FEM framework is to treat everything as a well-defined abstraction, simplifying the main assembly loop by delegating complex calculations to specialized objects.

-----

### The Core Philosophy of a Mature FEM Framework

A mature FEM core treats everything as a well-defined abstraction. The main assembly loop should be as simple as possible, delegating all complex calculations to specialized objects. The current `assemble` methods, with their internal logic for gathering DOFs and calculating Jacobians, are a sign that the `PhysicsField` class is doing too much work.

The goal is to move from this:
`PhysicsField` -> asks `DOFManager` for indices -> asks `ShapeFunctions` for values -> **calculates everything itself**

To this:
`PhysicsField` -> asks an `Element` for its pre-calculated `FEValues` -> **simply uses those values to fill the matrix**

-----

### **IV. Future Development Requirements**

**Congratulations!** You have successfully completed the core implementation of the vector field solver, including the necessary upgrades to `DOFManager` and `PhysicsField`, and have delivered a functional 3D magnetostatics solver with foundational tests. Our focus now shifts to in-depth optimization, theoretical expansion, and feature completeness.

### **1. Solver Performance Optimization & Robustness Enhancement**
* **Goal**: To significantly improve the solver's speed for large-scale problems, making it competitive with commercial software, and to ensure the absolute reliability of its results.
* **Key Tasks**:
    * **Advanced Preconditioning**: Research and implement more sophisticated preconditioners (e.g., Incomplete LU factorization `Eigen::IncompleteLUT`) for the `BiCGSTAB` solver to accelerate the convergence of complex electromagnetic problems.
    * **Parallelized Assembly**: Leverage the existing OpenMP integration to refactor the matrix assembly process (the element loop within the `assemble` function) for parallel execution and benchmark the performance gains.
    * **Comprehensive Benchmarking**: Establish a rigorous benchmark suite for the magnetostatics solver. Identify classic models with analytical solutions (e.g., Helmholtz coils) and write tests to perform **quantitative** accuracy validation, ensuring the error is within an acceptable tolerance.
    * **Enhanced Post-Processing**: Develop a new post-processor, similar to the `HeatFluxCalculator`, to compute the magnetic field **B** = ∇×**A**. This will be crucial for automated validation of the magnetic field results within the test suite.

### **2. Advancing into High-Frequency Applications: The Frequency-Domain Solver**
* **Goal**: To extend the solver's capabilities from static fields to time-harmonic fields (frequency domain), which is the core requirement for microwave, RF, and signal integrity simulations.
* **Governing Equation (Vector Helmholtz Equation)**: $\nabla \times \left( \frac{1}{\mu_r} \nabla \times \mathbf{E} \right) - k_0^2 \epsilon_r \mathbf{E} = 0$
    * Where **E** is the complex-valued electric field vector.
    * $k_0$ is the free-space wavenumber.
    * $\epsilon_r$ and $\mu_r$ are the relative permittivity and permeability.
* **Implementation Steps**:
    1.  **Support for Complex Arithmetic**: Extend the `LinearSolver` and `PhysicsField` classes to handle `std::complex<double>` matrices and vectors. The Eigen library provides native support for this.
    2.  **Create a `FrequencyDomainSolver`**: Develop a new `Solver` class designed to handle complex-valued linear systems.
    3.  **Implement `WavePropagation3D` Physics Field**: Create a new physics module to solve the Vector Helmholtz equation. This will require constructing complex-valued element stiffness matrices in the `assemble` method.
    4.  **Implement Advanced Boundary Conditions**: High-frequency simulations require **Absorbing Boundary Conditions (ABC)** or **Perfectly Matched Layers (PML)** to simulate open, non-reflecting boundaries. Begin by implementing a first- or second-order ABC.

### **3. Improving Project Generality and Usability**
* **Goal**: To enhance the project's modularity and ease of use, moving it closer to a general-purpose FEM solver core.
* **Key Tasks**:
    * **Non-Linear Material Support**: Extend the `Material` class to handle non-linear material properties, such as the B-H curve for magnetic materials. This will require implementing a **Newton-Raphson** iterative scheme within the solver.
    * **Advanced Meshing Capabilities**: Enhance the `Importer` to parse and utilize **Physical Groups** from Gmsh files. This will enable the application of different material properties and boundary conditions to distinct geometric regions within a single mesh.
    * **Documentation and Examples**: Create detailed documentation and clear example cases for newly implemented features, especially the frequency-domain solver and non-linear materials. This is critical for making the project understandable and usable by others.