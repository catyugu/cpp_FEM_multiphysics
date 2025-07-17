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

## **III. Reconstruction  Requirements**
Yes, the structure of the project is an **excellent foundation** and demonstrates a strong grasp of the principles behind a good FEM framework. However, to evolve into a truly **mature and scalable FEM core**, there are several key architectural shifts to consider.

Your current design excels at separating concerns at a high level—the `PhysicsField`, `Solver`, and `Core` namespaces are well-defined. The recent work on p-refinement and mixed-order solvers was a significant and necessary step forward.

The path from this strong foundation to a mature core involves pushing this abstraction principle even further down, primarily by making the **`Element` itself the central, intelligent entity** in the assembly process.

-----

### The Core Philosophy of a Mature FEM Framework

A mature FEM core treats everything as a well-defined abstraction. The main assembly loop should be as simple as possible, delegating all complex calculations to specialized objects. The current `assemble` methods, with their internal logic for gathering DOFs and calculating Jacobians, are a sign that the `PhysicsField` class is doing too much work.

The goal is to move from this:
`PhysicsField` -\> asks `DOFManager` for indices -\> asks `ShapeFunctions` for values -\> **calculates everything itself**

To this:
`PhysicsField` -\> asks an `Element` for its pre-calculated `FEValues` -\> **simply uses those values to fill the matrix**

-----

### The Reconstruction Guide: A 3-Phase Evolution

Here is a step-by-step guide to refactor your project into a more mature, extensible, and professional-grade architecture.

#### Phase 1: The "Smart Element" Refactoring (The Core Task)

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

#### Phase 2: Simplify the `assemble` Methods

With the new `FEValues` object, the `assemble` method in every single physics field becomes dramatically simpler, cleaner, and nearly identical.

Notice how the `PhysicsField` no longer knows anything about Jacobians or shape function derivatives. It just asks for the values it needs. **This is the hallmark of a mature FEM core.**

#### Phase 3: A Dedicated Post-Processing System

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


## **IV. Future Development Requirements**

With a robust p-refinement system in place, we can now focus on expanding the framework's capabilities. Please address the following tasks.

### **1. Implement Higher-Order (p > 2) Formulations**
* **Goal**: Extend the framework to support cubic, quartic, etc., elements.
* **Requirements**:
  * Implement the mathematical formulas for 3rd-order and higher shape functions in `utils/ShapeFunctions.cpp` for all element types.
  * Extend the `DOFManager` to handle DOFs on element **faces** (for 3rd-order 3D elements) and **internal/volume** DOFs (for 4th-order+ elements). This will likely require adding a `getFaceEquationIndex` method.
  * Update the `get_element_dofs` lambda function inside the `assemble` methods to correctly gather these new face and volume DOFs.

### **2. Implement Electromagnetic Field (Maxwell's Equations)**
* **Goal**: Add a full electromagnetic field simulation. The key challenge is upgrading the core components to handle vector-valued variables. A 3D Magnetostatics formulation is a great first step.
* **Guiding Equation**: $\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) = \mathbf{J}$
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

### **3. Add Support For Variable Material Properties**
* **Goal**: Simulate materials with variable properties (e.g., thermal conductivity changing with temperature).
* **Requirements**:
  * Add a `Material::getProperty` method that accepts a `time` argument.
  * Update the `PhysicsField::assemble` method to use the `time` argument to retrieve the current field values.
  * Create a new test case for a simple variable-property problem to validate the implementation.
  * If you have better constructions for this, please share!