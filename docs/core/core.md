# **Core Namespace**

This namespace contains the fundamental, high-level components that orchestrate the simulation framework. These classes manage the overall problem setup, degrees of freedom, material properties, boundary conditions, and domain sources.

---
## **Classes**

### **Problem**

* **Description**: The main orchestrator class. It owns all the components of a simulation (mesh, fields, solver) and controls the high-level workflow.
* **Public Functions**:
  * `Problem(std::unique_ptr<Mesh> mesh)`
  * `~Problem()`
  * `addField(std::unique_ptr<Physics::PhysicsField> field)`
  * `setup()`: This method now collects the element order for each physics field and passes it to the `DOFManager::build()` method, ensuring correct DOF mapping for higher-order elements.
  * `setIterativeSolverParameters(int max_iter, double tol)`: Sets parameters for iterative solvers.
  * `setTimeSteppingParameters(double time_step, double total_time)`: Sets parameters for transient simulations.
  * `setLinearSolverType(Solver::SolverType type)`: **New method** to specify the linear solver to be used (LU or BiCGSTAB).
  * `solveSteadyState()`
  * `solveTransient()`
  * `exportResults(const std::string& filename) const`
  * `getMesh() const`: returns `const Mesh&`
  * `getDofManager() const`: returns `const DOFManager&`
  * `getField(const std::string& var_name) const`: returns `Physics::PhysicsField*`
  * `getFields() const`: returns `const std::vector<std::unique_ptr<Physics::PhysicsField>>&`
  * `getCouplingManager()`: returns `CouplingManager&`
  * `getMaxIterations() const`, `getConvergenceTolerance() const`, `getTimeStep() const`, `getTotalTime() const`: Accessors for solver and time-stepping parameters.
  * `getLinearSolverType() const`: **New accessor** to retrieve the currently set linear solver type.
* **Private Members**:
  * `mesh_`: `std::unique_ptr<Mesh>`
  * `dof_manager_`: `std::unique_ptr<DOFManager>`
  * `solver_`: `std::unique_ptr<Solver::Solver>`
  * `fields_`: `std::vector<std::unique_ptr<Physics::PhysicsField>>`
  * `coupling_manager_`: `CouplingManager`
  * `max_iterations_`, `convergence_tolerance_`, `time_step_`, `total_time_`: Parameters for solvers.
  * `linear_solver_type_`: `Solver::SolverType` (Default: `LU`).

### **DOFManager**

* **Description**: The Degree of Freedom (DOF) Manager. It intelligently maps physics variables to unique indices in the global system of equations. It is designed to handle both standard vertex-based DOFs and the additional DOFs required for higher-order element formulations (e.g., nodes on element edges).
* **Public Functions**:
  * `DOFManager(Mesh& mesh)`
  * `registerVariable(const std::string& var_name)`
  * `build(const std::map<std::string, int>& field_orders)`: This method now requires a map of variable names to their corresponding element orders (`field_orders`). It uses this information to correctly build the DOF map, including assigning DOFs to higher-order edge nodes if required by any field.
  * `getEquationIndex(int node_id, const std::string& var_name) const`: Gets the global equation index for a **vertex** DOF.
  * `getEdgeEquationIndex(std::vector<int> node_ids, const std::string& var_name) const`: Gets the global equation index for a **higher-order edge** DOF. The edge is uniquely identified by the sorted IDs of its two vertex nodes.
  * `getNumEquations() const`: returns `size_t`
  * `getVariableNames() const`: returns `const std::vector<std::string>&`
* **Private Members**:
  * `mesh_`: `Mesh&`
  * `variable_names_`: `std::vector<std::string>`
  * `vertex_dof_map_`: `std::map<VertexDofKey, int>` (Maps vertex IDs and variable indices to global equation indices).
  * `edge_dof_map_`: `std::map<EdgeDofKey, int>` (New map for higher-order DOFs on edges).
  * `num_equations_`: `size_t`

### **Material**

* **Description**: Encapsulates the physical properties of a material. This class is designed to handle both constant and temperature-dependent properties.
* **Public Functions**:
  * `Material(const std::string& name)`
  * `setProperty(const std::string& prop_name, double value)`: Sets a constant scalar property.
  * `setTempDependentProperty(const std::string& prop_name, const std::map<std::string, double>& params)`: Sets parameters for a temperature-dependent property model (e.g., for electrical conductivity).
  * `getProperty(const std::string& prop_name) const`: Returns a constant property value.
  * `getProperty(const std::string& prop_name, double temperature) const`: Returns a property value evaluated at a given `temperature`. This overload enables temperature-dependent material properties.
  * `getName() const`: returns `const std::string&`
* **Private Members**:
  * `name_`: `std::string`
  * `properties_`: `std::map<std::string, std::any>` (Uses `std::any` to store both `double` for constant properties and `std::map<std::string, double>` for temperature-dependent model parameters).

### **BoundaryCondition**

* **Description**: An abstract base class for applying conditions on the boundary of the domain. It modifies the system of equations.
* **Public Functions**:
  * `BoundaryCondition(std::string tag = "")`
  * `apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const`: Pure virtual function to apply the BC.
  * `getEquationIndex() const`: **New pure virtual function** to retrieve the global equation index this BC applies to.
  * `getTag() const`: Returns the BC's tag.
* **Derived Classes**:
  * `DirichletBC`: Has two constructors: one for standard vertex nodes (taking a `node_id` and `var_name`) and a second, more direct constructor for higher-order DOFs (taking an `equation_index` directly). It uses the direct elimination method for application.
  * `NeumannBC`, `CauchyBC`.

### **SourceTerm**

* **Description**: An abstract base class for applying sources within the domain (e.g., a heat source). It modifies the force vector `F`.
* **Public Functions**:
  * `SourceTerm(std::string tag)`
  * `apply(Eigen::MatrixXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name) const`: Pure virtual function to apply the source.
  * `getTag() const`: Returns the source's tag.
* **Derived Classes**: `VolumetricSource`.