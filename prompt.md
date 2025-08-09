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

### Task 1: Decouple Physics Fields from Problem Definition (Dependency Injection)

-   **Goal**:
  -   To eliminate the direct dependency of `PhysicsField` subclasses on the `Problem` class, making the physics field modules independent, testable, and reusable "plugins".
-   **Motivation**:
  -   Currently, the `PhysicsField::assemble` method retrieves material properties via `problem_->getMaterial()`, creating tight coupling. An independent physics field should not need to know which specific "Problem" is using it. This aligns with the principles of modularity and reusability emphasized in the book.
-   **Implementation Strategy**:
  1.  Create a `MaterialManager` class as the material provider.
  2.  Modify the signature of the `assemble` method in the `PhysicsField` base class to accept a reference to the material manager, for instance:
      ```cpp
      // PhysicsField.hpp
      virtual void assemble(const MaterialManager& materials) = 0;
      ```
  3.  In the `Problem::solve()` method, invoke `physics_field_->assemble(*this);`, passing itself as the material manager.
  4.  Within the subclasses of `PhysicsField` (e.g., `Magnetic3D`), obtain the material for an element via the passed `materials` reference, not through the `problem_` member variable.

### Task 2: Implement Geometric Information Caching for Coupled-Field Analysis

-   **Goal**:
  -   To resolve the performance bottleneck in nonlinear iterative or transient analyses caused by the repeated calculation of unchanging geometric information.
-   **Motivation**:
  -   In many analyses (like the nonlinear electro-thermal coupling), the mesh geometry is fixed. Consequently, quantities like the `[B]` matrix, `detJ`, and shape function values `N` at integration points are constant for each element. Recomputing these in every iteration is a major source of performance degradation.
-   **Implementation Strategy**:
  1.  **Transform `FEValues` into a persistent object**: Modify the `Element` class to **own** a `std::unique_ptr<FEValues> fe_values_` member.
  2.  **Implement lazy loading/caching**: Implement a `getFEValues()` method in the `Element` class. This method checks if `fe_values_` has already been created. On the first call, it creates the `FEValues` object and runs its `reinit` method to perform the expensive one-time geometric calculations, storing the results (`B_matrices_`, `detJ_x_weights_`, etc.) as member variables of `FEValues`. Subsequent calls simply return the pointer to the already-cached `fe_values_` object.
  3.  **Refactor the coupling update function**: Modify functions like `ElectroThermalCoupling::updateSourceTerm` to no longer create a new `FEValues` instance in each iteration. Instead, they will retrieve the persistent, pre-computed object via `elem_ptr->getFEValues()`, drastically improving performance.

### Task 3: Implement the Element-by-Element (EBE) Iterative Solver

-   **Goal**:
  -   To implement a solver that does not rely on the explicit assembly of the global sparse matrix, breaking through the "memory wall" for large-scale problems.
-   **Motivation**:
  -   This is the core strategy for large problems presented in "Programming the Finite Element Method". Explicitly storing the global matrix becomes infeasible for systems with millions of DOFs. The EBE + Preconditioned Conjugate Gradient (PCG) method is the gold standard for solving this issue and is the foundation for high-performance parallel computing.
-   **Implementation Strategy**:
  1.  Create a new solver class, e.g., `EBE_PCG_Solver`.
  2.  Implement its core method, `matVecProduct(const Mesh& mesh, const Vector& p)`. This function will not access any global matrix. Instead, it will compute the global matrix-vector product by iterating through all elements, performing local `ke * pe` calculations, and accumulating the results into a global vector via a "Gather-Scatter" procedure.
  3.  Implement the main loop of the PCG algorithm, calling `matVecProduct` whenever the operation `A*p` is required.
  4.  Add an `assembleDiagonal()` method to the `LinearSystem` class to provide the simplest diagonal preconditioner for the PCG solver.
  5.  Integrate the new `EBESolver` into the `Problem` class as an optional solution strategy.
