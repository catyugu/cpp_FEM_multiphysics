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
  * `setup()`
  * `setIterativeSolverParameters(int max_iter, double tol)`
  * `setTimeSteppingParameters(double time_step, double total_time)`
  * `solveSteadyState()`
  * `solveTransient()`
  * `exportResults(const std::string& filename) const`
  * `getMesh() const`: returns `const Mesh&`
  * `getDofManager() const`: returns `const DOFManager&`
  * `getField(const std::string& var_name) const`: returns `Physics::PhysicsField*`
  * `getFields() const`: returns `const std::vector<std::unique_ptr<Physics::PhysicsField>>&`
  * `getCouplingManager()`: returns `CouplingManager&`
* **Private Members**:
  * `mesh_`: `std::unique_ptr<Mesh>`
  * `dof_manager_`: `std::unique_ptr<DOFManager>`
  * `solver_`: `std::unique_ptr<Solver::Solver>`
  * `fields_`: `std::vector<std::unique_ptr<Physics::PhysicsField>>`
  * `coupling_manager_`: `CouplingManager`

### **DOFManager**

* **Description**: The Degree of Freedom Manager. It maps the physics variables at each node to a unique index in the global system of equations.
* **Public Functions**:
  * `DOFManager(Mesh& mesh)`
  * `registerVariable(const std::string& var_name)`
  * `build()`
  * `getEquationIndex(int node_id, const std::string& var_name) const`
  * `getNumEquations() const`: returns `size_t`
  * `getVariableNames() const`: returns `const std::vector<std::string>&`

### **Material**

* **Description**: Encapsulates the physical properties of a material. This class is designed to handle both constant and temperature-dependent properties.
* **Public Functions**:
  * `Material(const std::string& name)`
  * `setProperty(const std::string& prop_name, double value)`
  * `setTempDependentProperty(const std::string& prop_name, const std::map<std::string, double>& params)`
  * `getProperty(const std::string& prop_name) const`: returns `double`
  * `getProperty(const std::string& prop_name, double temperature) const`: returns `double`
  * `getName() const`: returns `const std::string&`

### **BoundaryCondition**

* **Description**: An abstract base class for applying conditions on the boundary of the domain. It modifies the system of equations.
* **Public Functions**:
  * `BoundaryCondition(std::string tag = "")`
  * `apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const`: Pure virtual function to apply the BC.
  * `getTag() const`: Returns the BC's tag.
* **Derived Classes**: `DirichletBC`, `NeumannBC`, `CauchyBC`.

### **SourceTerm**

* **Description**: A new abstract base class for applying sources within the domain (e.g., a heat source). It modifies the force vector `F`.
* **Public Functions**:
  * `SourceTerm(std::string tag)`
  * `apply(Eigen::MatrixXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name) const`: Pure virtual function to apply the source.
  * `getTag() const`: Returns the source's tag.
* **Derived Classes**: `VolumetricSource`.