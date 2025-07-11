# **Solver Namespace**

This namespace is responsible for orchestrating the solution process of the simulation. It contains different solution strategies for single-field and coupled multiphysics problems.

---
## **Classes**

### **Solver**
* **Description**: This is an **interface** that defines the common structure for all solver implementations. It ensures that any solver can be used interchangeably by the main `Problem` class.
* **Public Functions**:
    * `~Solver()`: Virtual destructor.
    * `solveSteadyState(Core::Problem& problem)`: Pure virtual function to solve a steady-state problem.
    * `solveTransient(Core::Problem& problem)`: Pure virtual function to solve a transient (time-dependent) problem.

### **SingleFieldSolver**
* **Description**: A solver for problems involving only a single, uncoupled physical field (e.g., only heat transfer or only electrostatics).
* **Public Functions**:
    * `solveSteadyState(Core::Problem& problem) override`: Solves the steady-state system of equations for each physics field independently.
    * `solveTransient(Core::Problem& problem) override`: Solves the transient system of equations for a single physics field using a time-stepping scheme.

### **CoupledElectroThermalSolver**
* **Description**: A solver specifically designed for bidirectionally coupled electro-thermal problems, where electrical properties depend on temperature and Joule heating from the electrical field acts as a source for the thermal field.
* **Public Functions**:
    * `solveSteadyState(Core::Problem& problem) override`: Solves the steady-state coupled problem. It first solves the electric field, calculates the resulting Joule heat, and then uses that as a source to solve the heat transfer field.
    * `solveTransient(Core::Problem& problem) override`: Solves the time-dependent coupled problem. At each time step, it iteratively solves the electrical and thermal fields, updating material properties and heat sources until convergence is reached.

### **SolverFactory**
* **Description**: A factory class that creates and returns the appropriate solver object based on the physics fields present in the `Problem`.
* **Public Functions**:
    * `createSolver(const Core::Problem& problem)`: A static method that inspects the problem's fields. It returns a `CoupledElectroThermalSolver` if both "Voltage" and "Temperature" fields are present; otherwise, it returns a `SingleFieldSolver`.