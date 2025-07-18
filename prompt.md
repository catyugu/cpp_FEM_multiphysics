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

### The Reconstruction Guide: A 3-Phase Evolution

Here is a step-by-step guide to refactor your project into a more mature, extensible, and professional-grade architecture.

#### Phase 1: The "Smart Element" Refactoring (finished)

The goal is to make the `Element` class the primary source of all finite element calculation data. This removes the complex assembly logic from the individual physics fields.

**1. Create the `ElementGeometry` Class**
This class holds the pure geometric information of an element.

* **File**: `include/core/mesh/ElementGeometry.hpp`
* **Content**: It would store the raw vertex node coordinates and methods to compute things like edge lengths, face areas, and element volume. It knows nothing about shape functions or orders.

**2. Introduce `FEValues` (Finite Element Values)**
This is the most important new concept. The `FEValues` object is a "calculator" that, for a given quadrature point on an element, provides everything needed for assembly.

When an `FEValues` object is created, it pre-calculates and caches all shape function information for all quadrature points, making the assembly loop incredibly fast.

**3. Redesign the `Element` Class**
The `Element` class now becomes a container that holds its geometry and can create `FEValues` objects on demand.

* It would hold a `std::unique_ptr<ElementGeometry>`.
* Its primary new method would be: `std::unique_ptr<FEValues> create_fe_values(int order, int quad_order);`

#### Phase 2: Simplify the `assemble` Methods (finished)

With the new `FEValues` object, the `assemble` method in every single physics field becomes dramatically simpler, cleaner, and nearly identical.

Notice how the `PhysicsField` no longer knows anything about Jacobians or shape function derivatives. It just asks for the values it needs. **This is the hallmark of a mature FEM core.**

#### Phase 3: A Dedicated Post-Processing System (to do)

Mature software needs a flexible way to analyze results beyond just looking at the raw nodal values.

**1. Create a `PostProcessor` Base Class**

* **File**: `include/post/PostProcessor.hpp`
* **Key Method**: `virtual void compute_derived_quantities(const Problem& problem, const Eigen::MatrixXd& solution) = 0;`

**2. Implement Concrete Post-Processors**

* `HeatFluxCalculator`: Would take the temperature solution, compute the heat flux vector (**q** = -k∇T) at the quadrature points of each element, and add it as a new vector field to the results.
* `StrainCalculator`: For a solid mechanics problem, would compute the strain tensor from the displacement solution.

**3. Integrate with the `Problem` Class**

* Add a method `problem->addPostProcessor(...)`.
* After the solver finishes, the `Problem` class would loop through all registered post-processors and execute them, enriching the dataset before exporting the final `.vtk` file.

---

## **IV. Future Development Requirements**

With a robust p-refinement system in place, we can now focus on expanding the framework's capabilities. Please address the following tasks.

### **1. Add Support for Non-Linear Materials**
* **Goal**: Simulate materials whose properties change with the field itself (e.g., magnetic permeability changing with B-field strength).
* **Requirements**:
  * Document all public and private methods, classes, and functions in corresponding files in `docs/`.
  * Update the `Material::getProperty` method to accept the current element's field values as an optional argument.
  * Implement a non-linear solver loop (e.g., Newton-Raphson) within the `SingleFieldSolver` and `CoupledElectroThermalSolver`. This involves calculating a tangent stiffness matrix at each iteration.
  * Create a new test case for a simple non-linear problem to validate the implementation.


### **2. Implement Electromagnetic Field (Maxwell's Equations)**
* **Goal**: Add a full electromagnetic field simulation, initially focusing on a 3D Magnetostatics formulation as a key first step. The main challenge is upgrading the core components to handle **vector-valued variables**.
* **Guiding Equation (Magnetostatics)**: $\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) = \mathbf{J}$
    * **A** is the magnetic vector potential (a vector with components Ax, Ay, Az).
    * **μ** is the magnetic permeability.
    * **J** is the current density source.

* **Step-by-Step Implementation Guide**:
    1.  **Upgrade Core to Support Vector Fields**:
        * **`DOFManager`**:
            * Modify `registerVariable` to accept a component count: `void registerVariable(const std::string& var_name, int num_components = 1);`
            * Update the `build()` method to increment the equation counter by the number of components for each variable.
            * Ensure `getEquationIndex` returns the *starting* index for a variable's components (e.g., the index of `Ax`).
        * **`PhysicsField`**:
            * Add a new virtual method to the base class: `virtual int getNumComponents() const { return 1; }`. This defaults to 1 for all existing scalar fields.

    2.  **Implement the New `Magnetostatics3D` Physics Field**:
        * Create the new class inheriting from `Physics::PhysicsField`.
        * Override `getVariableName()` to return `"MagneticVectorPotential"`.
        * **Crucially, override `getNumComponents()` to return `3`**.
        * In the `assemble()` method, the local stiffness matrix `ke_local` for a 10-node quadratic tetrahedron will now be `30x30`. The "B-Matrix" must be reformulated to represent the **curl operator (`∇×`)** applied to vector-valued shape functions.

    3.  **Define a Current Source (`J`)**:
        * Create a new `PrescribedCurrentDensity` class that derives from `SourceTerm`.
        * Its constructor should accept a 3-component `Eigen::Vector3d` representing the current density.
        * Its `apply()` method must add the source contributions to the correct components of the global force vector `F`.

    4.  **Create a Validation Test**:
        * A great test case is a **long solenoid**.
        * Apply a circular `PrescribedCurrentDensity` in a cylindrical mesh to simulate the coil.
        * Set **A**=0 on the outer simulation boundary.
        * Validate that the computed magnetic field (**B** = ∇×**A**) is uniform inside the solenoid and matches the analytical solution.

### **3. Support for Advanced Research Problem Types (Solvers)**
* **Goal**: Expand the solver capabilities to handle more diverse research problems.
* **Requirements**:
    * **Time Domain Steady-State**: Implement a solver (or extend existing ones) to handle time-domain (transient) problems.
    * **Time Domain Transient**: Implement a solver (or extend existing ones) to handle time-domain (transient) problems.
    * **Frequency Domain Steady-State**: Implement a solver (or extend existing ones) to handle time-harmonic (AC) problems, where the solution is a complex number representing amplitude and phase. This will involve working with complex Eigen matrices.
    * **Coupled Transient Solvers**: Further enhance the robustness and efficiency of existing transient solvers, especially for highly non-linear and strongly coupled multiphysics problems.

### **4. Generalize Variable Material Properties for Coupled Fields**
* **Goal**: Extend the material property system to allow properties to depend on *any* relevant coupled field (not just temperature).
* **Current Status**: The `Material` class already supports temperature-dependent properties (e.g., electrical conductivity depending on temperature) via its `getProperty(const std::string& prop_name, double temperature) const` overload. This is already utilized in the `ElectroThermalCoupling`.
* **Requirements**:
    * **Flexible Property Evaluation**: Investigate a more generic mechanism within `Material::getProperty` to evaluate properties based on an arbitrary set of input field values (e.g., a map `std::map<std::string, double> field_values`). This would allow properties to depend on pressure, concentration, or other coupled variables.
    * **Physics Field Integration**: Update relevant `PhysicsField::assemble` methods to retrieve the necessary field values from the current solution of *other* fields (if coupled) and pass them to the material's property evaluation function.
    * **Testing**: Create new test cases that demonstrate material properties depending on various coupled fields.
