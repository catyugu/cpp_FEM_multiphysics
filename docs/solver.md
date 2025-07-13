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
* **Workflow (Steady-State)**:
  1.  Solves the electric field system.
  2.  Delegates to a `Coupling` object (e.g., `ElectroThermalCoupling`) to calculate the Joule heat and create the appropriate `SourceTerm` objects.
  3.  Solves the thermal field, which now includes the heat sources from the coupling step.
  4.  Repeats this process until the temperature solution converges.
  5.  It is responsible for stabilizing the system matrix of each field for the degrees of freedom of the *other* field during the solve.
* **Workflow (Transient)**:
  1.  **Time-Stepping Loop**: Advances the simulation from t=0 to `total_time` in increments of `time_step`.
  2.  **Inner Coupling Iteration**: At each time step, an inner iterative loop is performed to converge the coupled physics.
    * **Solve EMag Field**: The electrical field is solved based on the current (or previous inner iteration's) temperature solution.
    * **Execute Coupling**: The coupling mechanism calculates Joule heating based on the updated electrical solution.
    * **Solve Heat Field**: The thermal field is solved using the Backward Euler method: `(M/dt + K) * U_n = F + (M/dt) * U_{n-1}`. The force vector `F` includes the newly calculated Joule heating.
    * **Stabilization**: During both the EMag and Heat solves within the inner loop, the system matrix is stabilized by enforcing the degrees of freedom of the *other* physics field as known values.
    * **Inner Loop Convergence**: The inner loop continues until the relative change in the temperature solution falls below the `convergence_tolerance`.
  3.  **Update Previous Solution**: After the inner loop converges for a time step, the current temperature solution (`U_`) is saved as the previous solution (`U_prev_`) for the next overall time step.

### **SolverFactory**
* **Description**: A factory class that creates the appropriate solver (`SingleFieldSolver` or `CoupledElectroThermalSolver`) based on the physics fields present in the `Problem`.