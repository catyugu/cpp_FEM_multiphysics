# **Core Namespace**

This namespace contains the fundamental, high-level components that orchestrate the simulation framework. These classes manage the overall problem setup, degrees of freedom, material properties, and boundary conditions.

## **Classes**

### **Problem**

* **Description**: The main orchestrator class. It owns all the components of a simulation (mesh, fields, solver) and controls the high-level workflow.
* **Public Functions**:
  * Problem(std::unique\_ptr\<Mesh\> mesh)
  * \~Problem()
  * addField(std::unique\_ptr\<Physics::PhysicsField\> field)
  * setup()
  * setIterativeSolverParameters(int max\_iter, double tol)
  * setTimeSteppingParameters(double time\_step, double total\_time)
  * solveSteadyState()
  * solveTransient()
  * exportResults(const std::string& filename) const
  * getMesh() const: returns const Mesh&
  * getDofManager() const: returns const DOFManager&
  * getField(const std::string& var\_name) const: returns Physics::PhysicsField\*
  * getFields() const: returns const std::vector\<std::unique\_ptr\<Physics::PhysicsField\>\>&
  * getCouplingManager(): returns CouplingManager&
* **Private Members**:
  * mesh\_: std::unique\_ptr\<Mesh\>
  * dof\_manager\_: std::unique\_ptr\<DOFManager\>
  * solver\_: std::unique\_ptr\<Solver::Solver\>
  * fields\_: std::vector\<std::unique\_ptr\<Physics::PhysicsField\>\>
  * coupling\_manager\_: CouplingManager

### **DOFManager**

* **Description**: The Degree of Freedom Manager. It maps the physics variables at each node to a unique index in the global system of equations.
* **Public Functions**:
  * DOFManager(Mesh& mesh)
  * registerVariable(const std::string& var\_name)
  * build()
  * getEquationIndex(int node\_id, const std::string& var\_name) const
  * getNumEquations() const: returns size\_t
  * getVariableNames() const: returns const std::vector\<std::string\>&
* **Private Members**:
  * mesh\_: Mesh&
  * variable\_names\_: std::vector\<std::string\>
  * dof\_map\_: std::map\<std::pair\<int, int\>, int\>
  * num\_equations\_: size\_t

### **Material**

* **Description**: Encapsulates the physical properties of a material. This class is designed to handle both constant and temperature-dependent properties.
* **Public Functions**:
  * Material(const std::string& name)
  * setProperty(const std::string& prop\_name, double value)
  * setTempDependentProperty(const std::string& prop\_name, const std::map\<std::string, double\>& params)
  * getProperty(const std::string& prop\_name) const: returns double
  * getProperty(const std::string& prop\_name, double temperature) const: returns double
  * getName() const: returns const std::string&
* **Private Members**:
  * name\_: std::string
  * properties\_: std::map\<std::string, std::any\>

### **BoundaryCondition**

* **Description**: A polymorphic system for applying boundary conditions. This is an abstract base class.
* **Public Functions**:
  * apply(Eigen::SparseMatrix\<double\>& K, Eigen::MatrixXd& F) const: Apply the boundary condition to the global system of equations.
* **Derived Classes**:
  * DirichletBC
  * NeumannBC
  * CauchyBC

### **DirichletBC**

* **Description**: Fixes a degree of freedom to a specific value.
* **Public Functions**:
  * DirichletBC(const DOFManager& dof\_manager, int node\_id, const std::string& var\_name, Eigen::VectorXd value)
  * apply(Eigen::SparseMatrix\<double\>& K, Eigen::MatrixXd& F) const override

### **NeumannBC**

* **Description**: Applies a known flux to a node, modifying the force vector F.
* **Public Functions**:
  * NeumannBC(const DOFManager& dof\_manager, int node\_id, const std::string& var\_name, Eigen::VectorXd flux\_value)
  * apply(Eigen::SparseMatrix\<double\>& K, Eigen::MatrixXd& F) const override

### **CauchyBC**

* **Description**: Models mixed/convective boundaries, modifying both K and F.
* **Public Functions**:
  * CauchyBC(const DOFManager& dof\_manager, int node\_id, const std::string& var\_name, Eigen::VectorXd h, Eigen::VectorXd T\_inf)
  * apply(Eigen::SparseMatrix\<double\>& K, Eigen::MatrixXd& F) const override