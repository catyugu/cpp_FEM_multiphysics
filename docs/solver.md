# **Solver Namespace**

This namespace is responsible for orchestrating the solution process of the simulation. It contains different solution strategies for single-field and coupled multiphysics problems.

---
## **Classes**

### **Solver**
* **Description**: An **interface** that defines the common structure for all solver implementations.

### **LinearSolver**
* **Description**: A static utility class that solves a linear system of equations using the Eigen library. It provides access to both direct (`SparseLU`) and iterative (`BiCGSTAB`) solvers.

### **SingleFieldSolver**
* **Description**: A solver for problems involving only a single, uncoupled physical field.

### **CoupledElectroThermalSolver**
* **Description**: A solver for bidirectionally coupled electro-thermal problems. It uses an iterative approach to find a converged solution between the electric and thermal fields.
* **Key Features**:
  * **Mixed-Order Capability**: The solver's stabilization logic is now robust enough to correctly handle coupled simulations where the electric and thermal fields are discretized with different mathematical orders (e.g., coupling a 2nd-order Voltage field with a 1st-order Temperature field).
  * **Transient and Steady-State**: Implements the iterative workflow for both steady-state and transient problems, ensuring stability at each step.
  * **Coupling Mechanism**: Delegates to `Coupling` objects to manage the physics interactions, such as calculating Joule heat.

### **SolverFactory**
* **Description**: A factory class that creates the appropriate solver (`SingleFieldSolver` or `CoupledElectroThermalSolver`) based on the physics fields present in the `Problem`.