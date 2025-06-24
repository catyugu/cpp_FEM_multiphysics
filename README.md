# C++ FEM Multiphysics
## Description
This project is a C++ implementation of a finite element method (FEM) for solving a multiphysics problem.
## Requirements
* C++17
* CMake >= 3.17
* Eigen >= 3.3
* GTest >= 1.10
## Build Guide
```bash
mkdir cmake-build-debug && cd cmake-build-debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
./run_tests ## Run google tests
./cpp_FEM_multiphysics ## Run main program
```
## Current Functions
* Finite Element Method Core.
* Solving 1D Heat&ElectroMagnetic Field and their coupling.
* Solving 2D Heat&ElectroMagnetic Field and their coupling.
* Export result in vtk file.
## Core File Structure
```bash
cpp_FEM_multiphysics
├── CMakeLists.txt      # Main build script
├── include/            # Header files (.hpp)
│   ├── core/           # Core architectural components
│   ├── io/             # Input/Output utilities (e.g., Exporter)
│   ├── physics/        # Physics-specific modules
│   └── utils/          # General utilities (e.g., Logger)
├── src/                # Source files (.cpp)
│   ├── core/
│   ├── io/
│   ├── physics/
│   └── main.cpp        # Main application entry point
└── tests/              # Google Test source files
```
Of course. Here is the content of the `README.md` document in plain text.

## Architecture and Class Reference

The solver is built on an object-oriented design that separates concerns into distinct architectural layers: Core, Physics, and IO.

### `core` Namespace

This namespace contains the fundamental building blocks of the FEM framework.

#### `Problem`

* **Description:** The main orchestrator class. It owns all the components of a simulation (mesh, fields, BCs) and controls the high-level workflow.
* **Key Methods:**
    * `Problem(std::unique_ptr<Mesh> mesh)`: Constructor that takes ownership of the mesh.
    * `addField(std::unique_ptr<Physics::PhysicsField> field)`: Adds a physics simulation (e.g., Heat2D) to the problem.
    * `setSolverParameters(int max_iter, double tol)`: Configures the iterative solver.
    * `setup()`: Initializes all components, builds the DOF map.
    * `solve()`: Runs the simulation. Contains the logic for uncoupled, weakly-coupled, and strongly-coupled iterative solves.
    * `exportResults(const std::string& filename)`: Exports the final results to a file.

#### `Mesh`, `Node`, `Element`

* **Description:** These classes define the problem geometry.
* **`Node`:** Represents a point in space with an ID and coordinates.
* **`Element`:** An abstract base class for all mesh elements.
* **`LineElement`, `TriElement`:** Concrete implementations for 1D and 2D elements. They contain geometric calculations like length, area, and the B-matrix (strain-displacement matrix).
* **`Mesh`:** A container for all nodes and elements in the simulation. Includes static factory methods (`create_uniform_1d_mesh`, `create_uniform_2d_mesh`) to generate simple meshes.

#### `Material`

* **Description:** Encapsulates the physical properties of a material. This class is designed to handle both constant and temperature-dependent properties.
* **Key Methods:**
    * `setProperty(name, value)`: Sets a simple, constant property.
    * `setTempDependentProperty(name, params)`: Sets the parameters for a functional material model (e.g., a linear model for conductivity).
    * `getProperty(name)`: Gets a constant property value.
    * `getProperty(name, temperature)`: Gets a property value at a specific temperature, evaluating the underlying model.

#### `BoundaryCondition`

* **Description:** A polymorphic system for applying boundary conditions.
* **`BoundaryCondition`:** The abstract base class. Its key method is `apply(K, F)`, which modifies the global stiffness matrix and force vector.
* **`DirichletBC`:** Fixes a degree of freedom to a specific value.
* **`NeumannBC`:** Applies a known flux to a node, modifying the force vector `F`.
* **`CauchyBC`:** Models mixed/convective boundaries, modifying both `K` and `F`.

#### `DOFManager`

* **Description:** The Degree of Freedom Manager. It maps the physics variables at each node (e.g., "Temperature at Node 5") to a unique index (row/column) in the global system of equations.

#### `LinearSolver`

* **Description:** A simple, static wrapper class around the Eigen library's sparse solver (`SparseLU`). It provides a single function `solve(A, b, x)`.

### `physics` Namespace

This namespace contains the specific implementations for different physical simulations.

#### `PhysicsField`

* **Description:** An abstract base class for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module.
* **Key Methods:**
    * `setup(...)`: Initializes the field, resizing matrices and vectors.
    * `assemble()`: Contains the core FEM logic to compute the element stiffness matrices and assemble them into the global system.
    * `getSolution()`: Returns the vector of solved nodal values.
    * `addBC(...)`: Attaches a boundary condition to the field.

#### Concrete Physics Classes (`Heat1D`, `Heat2D`, `EMag1D`, `EMag2D`)

* **Description:** These classes implement the `assemble` method for their specific physics. For example, `Heat2D` implements the B-matrix formulation (`Area * B^T * D * B`) for 2D triangular elements. Coupled classes also include methods to calculate or receive data from other fields (e.g., `EMag2D::calculateJouleHeat`, `Heat2D::setVolumetricHeatSource`).

### `io` Namespace

#### `Exporter`

* **Description:** Handles the writing of simulation results to files.
* **Key Methods:**
    * `write_vtk(filename, problem)`: A static method that takes a `Problem` object and writes its mesh and all solved nodal data into a single, standard `.vtk` file for visualization.