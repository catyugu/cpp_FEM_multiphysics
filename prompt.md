# **C++ FEM Multiphysics Project: Contributor Onboarding Guide (Phase 2\)**

Welcome to the next phase of the **C++ FEM Multiphysics** project\! This document serves as your guide to understanding the existing architecture, adhering to our coding practices, and tackling the next set of features. Our goal is to continue building a robust, scalable, and maintainable finite element analysis framework.

## **I. Onboarding and Project Philosophy**

This project is a C++ implementation of the finite element method (FEM) designed to solve coupled multiphysics problems. The core philosophy is to create a modular system where new physics, element types, and solution strategies can be added with minimal changes to the existing framework.

We prioritize:

* **Modularity**: Separating core logic from physics-specific implementations.
* **Clarity**: Writing self-documenting code with clear and consistent naming.
* **Correctness**: Ensuring all implementations are backed by rigorous testing against analytical solutions or benchmarks.
* **Safety**: Implementing robust error handling and exception safety.
* **Extensibility**: Allowing new physics, element types, and solution strategies to be added easily.
* **Testing**: All new features must be accompanied by comprehensive tests.
* **Documentation**: All new code must be documented in the corresponding file in the docs/ folder.

## **II. Coding Guidelines and Practices**

Adherence to these guidelines is essential for maintaining code quality and consistency.

### **1\. Project Structure**

The repository is organized into distinct directories. Please place new files in their appropriate locations.

cpp\_FEM\_multiphysics/  
├── CMakeLists.txt      \# Main build script  
├── docs/               \# High-level documentation  
├── include/            \# Header files (.hpp)  
│   ├── core/           \# Core architectural components (Problem, Mesh, Solver)  
│   ├── io/             \# Importer/Exporter utilities  
│   ├── physics/        \# Physics-specific modules (Heat2D, Current1D)  
│   └── utils/          \# General utilities (SimpleLogger, Exceptions)  
├── src/                \# Source files (.cpp)  
│   ├── core/  
│   ├── io/  
│   ├── physics/  
│   └── main.cpp        \# Main application entry point  
└── tests/              \# Google Test source files

### **2\. C++ Language and Style**

* **C++ Standard**: The project uses **C++17**. Use C++17 features where they improve clarity and safety.
* **Dependencies**: The primary dependencies are **Eigen** for linear algebra and **GTest** for testing.
* **Naming Conventions**:
    * **Namespaces**: PascalCase (e.g., Core, Physics, IO). All code should reside within a namespace.
    * **Classes/Structs**: PascalCase (e.g., DOFManager, TriElement).
    * **Functions/Methods**: camelCase (e.g., getEquationIndex, solveSteadyState).
    * **Member Variables**: snake\_case\_ with a trailing underscore (e.g., mesh\_, num\_equations\_). This is a strict project rule.
    * **Local Variables**: snake\_case (e.g., k\_triplets, num\_steps).

### **3\. Core Architectural Principles**

* **Problem Class**: The Core::Problem class is the main orchestrator. It owns the mesh, DOF manager, and physics fields. It does **not** contain direct solving logic.
* **Solver Strategy**: The actual solution logic is delegated to Solver objects (e.g., SingleFieldSolver, CoupledElectroThermalSolver). The correct solver is chosen by the SolverFactory.
* **PhysicsField Abstraction**: All physics (e.g., thermal, electrical) must inherit from the Physics::PhysicsField base class. This class defines the interface for assembling matrices (K\_, M\_) and the RHS vector (F\_).
* **Memory Management**: The project uses std::unique\_ptr to manage the lifetime of major components owned by the Problem class. Raw pointers are used within the Mesh class for nodes and elements, which are treated as non-owning containers.

## **III. Future Development Requirements**

The following are high-priority areas for the next development cycle. Please address them in the order presented.

### **1\. Complete 3D Physics Capabilities**

* **Goal**: Extend the framework to fully support 3D simulations for electromagnetics, matching the existing 3D heat transfer capabilities.
* **Tasks**:
    1. **Implement Current3D**: Create a Physics::Current3D class that inherits from PhysicsField. The primary variable will be "Voltage". The assemble method should compute the element stiffness matrix for TetElement using its getBMatrix() and getVolume() methods.
    2. **Implement Magnetic3D**: Create a Physics::Magnetic3D class. The primary variable will be "MagneticPotential". The assembly will be similar to Current3D but will use magnetic permeability from the Material class.
    3. **Validation Tests**: For both Current3D and Magnetic3D, create new test files in the /tests directory. Each test should solve a simple problem with a known analytical solution on a 3D mesh (e.g., a cube with fixed boundary conditions) to validate the implementation.

### **2\. Code Health and Refactoring**

* **Goal**: Improve code maintainability and remove obsolete code.
* **Tasks**:
    1. **Review and Clean**: Perform a comprehensive review of all classes. Identify and remove any member variables, functions, or private helper methods that are no longer used or have been superseded by newer implementations.
    2. **Consolidate Redundancy**: Look for opportunities to consolidate redundant code, especially within the assemble methods of the various PhysicsField classes.

### **3\. Foundational Work for Higher-Order Elements**

* **Goal**: Prepare the framework to support more advanced, higher-order (e.g., quadratic) finite elements. This will significantly improve solution accuracy without requiring extremely fine meshes.
* **Tasks**:
    1. **Element Degree Control**:
        * Add a new integer member variable, element\_order\_, to the PhysicsField base class, with a corresponding setter (e.g., setElementOrder(int order)). The default value should be 1 (for linear elements).
    2. **Numerical Integration Utilities**:
        * Create a new utility class, Utils::Quadrature, to handle numerical integration.
        * This class should provide static methods to get integration points (Gauss points) and their corresponding weights for different element types (e.g., getLineQuadrature(int order), getTriangleQuadrature(int order)). Initially, implement this for linear elements (order \= 1).
    3. **Prepare for Higher-Order Shape Functions**:
        * The current Element classes (e.g., TriElement, TetElement) assume simple linear behavior. While you do not need to implement full quadratic elements yet, refactor the assemble methods in the physics classes to use a numerical integration loop. This loop will iterate over the quadrature points provided by the new Quadrature utility.
        * This change will replace the simple Area \* B^T \* D \* B calculation with a more general form: k\_e \= integral(B^T \* D \* B) dV, which is approximated by summing (B^T \* D \* B \* weight \* detJ) at each integration point. For now, with linear elements, the result will be the same, but the code structure will be ready for higher-order shape functions and their derivatives.