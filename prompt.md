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
│   │   ├── mesh/
│   │   └── sources/    \# SourceTerm classes
│   ├── io/             \# Importer/Exporter utilities
│   ├── physics/        \# Physics-specific modules (Heat2D, Current3D)
│   ├── solver/
│   ├── post/
│   └── utils/          \# General utilities (Logger, Quadrature)
├── src/                \# Source files (.cpp)
│   ├── core/
│   │   ├── bcs/
│   │   ├── coupling/
│   │   ├── mesh/
│   │   └── sources/
│   ├── io/             
│   ├── physics/    
│   ├── solver/
│   ├── post/
│   └── utils/
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

The framework has evolved into a sophisticated, modular system capable of handling complex multiphysics simulations with higher-order approximations. The current architecture centers around several key components that work together to provide a flexible and extensible platform for finite element analysis.

### **1. Core Component Architecture**

The framework is organized around these foundational components:

* **Problem**: Acts as the central coordinator, managing the physics fields, coupling mechanisms, and solution process.
* **Mesh**: Represents the discretized domain with support for various element types and refinement strategies.
* **Element**: Implements specific finite element formulations (triangular, tetrahedral, etc.) with support for different approximation orders.
* **Material**: Encapsulates physical properties that can be constant, spatially-varying, or dependent on solution fields.
* **PhysicsField**: Abstract base class for different physics modules (heat, current, magnetism) that defines the assembly of element matrices and vectors.

### **2. Advanced Degree of Freedom (DOF) Management**

To support p-refinement and higher-order approximations, the `DOFManager` has been completely redesigned:

* **Multiple DOF Maps**: It now maintains separate maps for different types of DOFs: one for vertices (`vertex_dof_map_`) and one for edges (`edge_dof_map_`).
* **Unique Edge Identification**: To uniquely identify a higher-order DOF on an edge, it uses a key composed of the **sorted** vertex IDs of that edge. This guarantees that `(Node1, Node2)` and `(Node2, Node1)` refer to the same edge DOF.
* **Order-Aware Build Process**: The `build()` method now requires the element orders of all physics fields to determine exactly which higher-order DOFs need to be created.

### **3. Sophisticated FEM Computation Engine**

The framework implements advanced FEM calculation capabilities:

* **FEValues**: A class that computes and caches shape functions, their derivatives, Jacobian determinants, and other geometric quantities at quadrature points.
* **Quadrature**: Provides integration rules of various orders for accurate numerical integration.
* **Shape Functions**: Supports hierarchical shape functions for arbitrary-order approximations, enabling p-refinement.

### **4. Robust Boundary Condition (BC) System**

Applying constraints to higher-order elements required a more sophisticated BC system:

* **Dual-Constructor `DirichletBC`**: The `DirichletBC` class now has two constructors: one that uses a `node_id` (for vertices) and a second that directly accepts an `equation_index`. This second constructor is essential for constraining the "virtual" higher-order DOFs, which do not have a node ID.
* **Consolidated Application**: The `PhysicsField::applyBCs()` method is now intelligent. It first consolidates all BCs into a map where each DOF index appears only once. This prevents common errors from redundantly applying constraints to the same node (e.g., a corner node shared by multiple boundary edges), which previously caused solver failures.
* **Tag-Based System**: BCs and sources can be applied using named tags, making the interface more intuitive and less error-prone.

### **5. Physics Module System**

The framework supports various physics modules, each implementing specific element formulations:

* **Heat2D/Heat3D**: For thermal conduction analysis with support for steady-state and transient solutions.
* **Current2D/Current3D**: For electric current flow analysis.
* **Magnetic3D**: For 3D magnetostatic analysis with vector field approximation.

### **6. Coupling Mechanisms**

A flexible system for handling multiphysics interactions:

* **Coupling**: Abstract base class defining the interface for field interactions.
* **ElectroThermalCoupling**: Implementation for coupled electric-thermal analysis with Joule heating.
* **SourceTerm**: Base class for representing various physical source terms in equations.

### **7. Solver and Performance**

The framework includes several solver strategies:

* **DirectSolver**: Uses sparse direct solvers for robust solution of linear systems.
* **Mixed-Order Coupled Solver**: The `CoupledElectroThermalSolver` has been upgraded to correctly stabilize and solve systems where different physics fields have different element orders.
* **Transient Solver**: Implements time-stepping schemes for dynamic analyses.

### **8. Input/Output System**

The framework provides tools for data import and export:

* **MeshImporter**: Reads mesh data from various formats (currently supporting MPHTXT).
* **ResultExporter**: Writes solution data to visualization formats (VTK/VTU).

This architecture provides a solid foundation for implementing complex multiphysics simulations while maintaining modularity, extensibility, and performance.

---

### **IV. Future Development Requirements**

**Congratulations!** You have successfully completed the core implementation of the vector field solver, including the necessary upgrades to `DOFManager` and `PhysicsField`, and have delivered a functional 3D magnetostatics solver with foundational tests. Our focus now shifts to in-depth optimization, theoretical expansion, and feature completeness.

### Task 1: Decouple Physics Fields from Problem Definition (Dependency Injection)

-   **Goal**:
  -   To eliminate the direct dependency of `PhysicsField` subclasses on the `Problem` class, making the physics field modules independent, testable, and reusable "plugins".
-   **Motivation**:
  -   Currently, the `PhysicsField::assemble` method retrieves material properties via `problem_->getMaterial()`, creating tight coupling. An independent physics field should not need to know which specific "Problem" is using it. This aligns with the principles of modularity and reusability emphasized in the book.
-   **Implementation Strategy**:
  1.  Create a `MaterialManager` class or use the `Problem` class itself as the material provider.
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
