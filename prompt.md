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

* **Mixed-Order Coupled Solver**: The `CoupledElectroThermalSolver` has been upgraded to correctly stabilize and solve systems where different physics fields have different element orders. This is crucial for advanced multiphysics simulations.

## **III. Future Development Requirements**

With a robust p-refinement system in place, we can now focus on expanding the framework's capabilities. Please address the following tasks.

### **1. Implement Higher-Order (p > 2) Formulations**
* **Goal**: Extend the framework to support cubic, quartic, etc., elements.
* **Requirements**:
  * Implement the mathematical formulas for 3rd-order and higher shape functions in `utils/ShapeFunctions.cpp` for all element types.
  * Extend the `DOFManager` to handle DOFs on element **faces** (for 3rd-order 3D elements) and **internal/volume** DOFs (for 4th-order+ elements). This will likely require adding a `getFaceEquationIndex` method.
  * Update the `get_element_dofs` lambda function inside the `assemble` methods to correctly gather these new face and volume DOFs.

### **2. Add Support for Non-Linear Materials**
* **Goal**: Simulate materials whose properties change with the field itself (e.g., magnetic permeability changing with B-field strength).
* **Requirements**:
  * Update the `Material::getProperty` method to accept the current element's field values as an optional argument.
  * Implement a non-linear solver loop (e.g., Newton-Raphson) within the `SingleFieldSolver` and `CoupledElectroThermalSolver`. This involves calculating a tangent stiffness matrix at each iteration.
  * Create a new test case for a simple non-linear problem to validate the implementation.

### **3. Performance Optimization with Parallelism**
* **Goal**: Significantly speed up the assembly process, which is the next major bottleneck after the linear solve.
* **Requirements**:
  * The `assemble` loop over elements is a prime candidate for parallelization. Use **OpenMP** (which is already included in the project) to parallelize this loop.
  * Ensure that the accumulation of triplets into the `k_triplets` and `m_triplets` vectors is done in a thread-safe manner. A standard approach is to have each thread work on its own private list of triplets and merge them all at the end.
  * Benchmark the performance improvement on a large test case.
