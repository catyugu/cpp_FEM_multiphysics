# **Solver Namespace**

This namespace is responsible for orchestrating the solution process of the simulation. It contains different solution strategies for single-field, coupled, steady-state, and transient problems.

---
## **Classes**

### **Solver**
* **Description**: An abstract base class that defines the interface for all solver types.
* **Public Functions**:
    * `solveSteadyState(Core::Problem& problem)`: Pure virtual function for steady-state analysis.
    * `solveTransient(Core::Problem& problem)`: Pure virtual function for time-dependent analysis.

### **SingleFieldSolver**
* **Description**: A solver for problems involving only a single, uncoupled physical field, or multiple fields that do not interact.
* **Key Logic**:
    * `solveSteadyState`: Iterates through each enabled physics field, calls its `assemble()`, `applySources()`, and `applyBCs()` methods, and then solves the linear system `K*U = F`.
    * `solveTransient`: For each time step, it iterates through each field and solves the time-discretized equation using a Backward Euler scheme: `(M/dt + K) * U_n+1 = F + (M/dt) * U_n`.

### **CoupledElectroThermalSolver**
* **Description**: A sophisticated solver for bidirectionally coupled electro-thermal problems, where electrical conductivity depends on temperature and Joule heating acts as a heat source.
* **Key Logic**:
    * **`solveSteadyState`**: Implements a robust, damped, Newton-like iterative scheme:
        1.  In each iteration, first solve the linear electrical problem using the temperature from the previous iteration to evaluate conductivity.
        2.  Call the `CouplingManager` to calculate the Joule heat source based on the new voltage solution.
        3.  Solve the heat transfer problem.
        4.  Check for convergence of the temperature field. A damping factor is used to improve stability.
        5.  A key stability enhancement is the "stabilization" of non-active DOFs. When solving for one field (e.g., voltage), the matrix rows/columns corresponding to the other field's (e.g., temperature) DOFs are temporarily modified to prevent a singular matrix, ensuring a robust solution even with strong coupling.
    * **`solveTransient`**: Implements a time-stepping loop. Within each time step, it performs an inner iteration loop similar to the steady-state solver to resolve the non-linear coupling at that specific point in time.

### **LinearSolver**
* **Description**: A static utility class that provides a unified interface to Eigen's linear solvers.
* **Public Functions**:
    * `solve(...)`: A static method that takes a sparse matrix `A`, a vector `b`, and solves `Ax = b`. It can be configured to use either a direct solver (`Eigen::SparseLU`) for smaller, robust solves, or an iterative solver (`Eigen::BiCGSTAB`) for larger systems where performance is critical.

### **SolverFactory**
* **Description**: A factory class that automatically creates the appropriate solver based on the physics fields present in the `Problem`.
* **Functionality**: It checks if both "Voltage" and "Temperature" fields are registered. If so, it returns a `std::unique_ptr<CoupledElectroThermalSolver>`; otherwise, it returns a `std::unique_ptr<SingleFieldSolver>`.