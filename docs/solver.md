# **Solver Namespace**

This namespace is responsible for orchestrating the solution process of the simulation. It contains different solution strategies for single-field and coupled multiphysics problems.

---
## **Classes**

### **Solver**
* **Description**: An **interface** that defines the common structure for all solver implementations.

### **LinearSolver**
* **Description**: A static utility class that solves a linear system of equations using the Eigen library. It provides access to both direct (`SparseLU`) and iterative (`BiCGSTAB`) solvers.
* **Key Function**:
  * `static void solve(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& b, Eigen::MatrixXd& x, SolverType solver_type = SolverType::LU, int max_iterations = 1000, double tolerance = 1e-9)`: Solves the linear system Ax = b. It now explicitly accepts `solver_type`, `max_iterations`, and `tolerance` parameters, which are passed from the `Problem` object.
* **Enum**:
  * `SolverType { LU, BiCGSTAB }`: Defines the available linear solver types.

### **SingleFieldSolver**
* **Description**: A solver for problems involving only a single, uncoupled physical field.
* **Key Update**: The `solveSteadyState` and `solveTransient` methods now correctly pass the linear solver type, maximum iterations, and convergence tolerance from the `Problem` object to the `LinearSolver::solve` method.

### **CoupledElectroThermalSolver**
* **Description**: A solver for bidirectionally coupled electro-thermal problems. It uses an iterative approach to find a converged solution between the electric and thermal fields.
* **Key Features**:
  * **Mixed-Order Capability**: The solver's stabilization logic is robust enough to correctly handle coupled simulations where the electric and thermal fields are discretized with different mathematical orders (e.g., coupling a 2nd-order Voltage field with a 1st-order Temperature field).
  * **Transient and Steady-State**: Implements the iterative workflow for both steady-state and transient problems, ensuring stability at each step.
  * **Coupling Mechanism**: Delegates to `Coupling` objects to manage the physics interactions, such as calculating Joule heat.

### **SolverFactory**
* **Description**: A factory class that creates the appropriate solver (`SingleFieldSolver` or `CoupledElectroThermalSolver`) based on the physics fields present in the `Problem`.