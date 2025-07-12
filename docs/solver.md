# **Solver Namespace**

This namespace is responsible for orchestrating the solution process of the simulation. It contains different solution strategies for single-field and coupled multiphysics problems.

---
## **Classes**

### **Solver**
* **Description**: An **interface** that defines the common structure for all solver implementations.

### **LinearSolver**
* **Description**: A static utility class that solves a linear system of equations using the Eigen library.

### **SingleFieldSolver**
* **Description**: A solver for problems involving only a single, uncoupled physical field. Its workflow is to call `assemble()`, `applySources()`, and `applyBCs()` for each field before solving.

### **CoupledElectroThermalSolver**
* **Description**: A solver for bidirectionally coupled electro-thermal problems. It uses an iterative approach to find a converged solution between the electric and thermal fields.
* **Workflow**:
  1.  Solves the electric field system.
  2.  Delegates to a `Coupling` object (e.g., `ElectroThermalCoupling`) to calculate the Joule heat and create the appropriate `SourceTerm` objects.
  3.  Solves the thermal field, which now includes the heat sources from the coupling step.
  4.  Repeats this process until the temperature solution converges.
  5.  It is responsible for stabilizing the system matrix of each field for the degrees of freedom of the *other* field during the solve.

### **SolverFactory**
* **Description**: A factory class that creates the appropriate solver (`SingleFieldSolver` or `CoupledElectroThermalSolver`) based on the physics fields present in the `Problem`.