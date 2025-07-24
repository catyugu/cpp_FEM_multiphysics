# **Solver Namespace**

This namespace is responsible for orchestrating the solution process of the simulation.

---
## **Classes**

### **SingleFieldSolver**
* **Description**: A solver for problems involving only a single, uncoupled physical field.
* **Key Update**: The `solveSteadyState` and `solveTransient` methods are now responsible for the order of operations. They first call `assemble()` on the field, then `applySources()`, and finally `applyBCs()`. This ensures that source terms are correctly added to the system before boundary conditions are imposed.

### **LinearSolver**
* **Description**: A static utility class that solves a linear system of equations using the Eigen library. It provides access to both direct (`SparseLU`) and iterative (`BiCGSTAB`) solvers.

### **CoupledElectroThermalSolver**
* **Description**: A solver for bidirectionally coupled electro-thermal problems.

### **SolverFactory**
* **Description**: A factory class that creates the appropriate solver based on the physics fields present in the `Problem`.
