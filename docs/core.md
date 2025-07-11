# **Core Namespace**

This namespace contains the fundamental building blocks of the FEM framework.

---
## **Classes**

### **BoundaryCondition**
* **Description**: A polymorphic system for applying boundary conditions. This is an abstract base class.
* **Public Functions**:
  * `apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const`: Apply the boundary condition to the global system of equations.
* **Derived Classes**:
  * `DirichletBC`
  * `NeumannBC`
  * `CauchyBC`

### **DirichletBC**
* **Description**: Fixes a degree of freedom to a specific value.
* **Public Functions**:
  * `DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value)`
  * `apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override`
* **Private Members**:
  * `equation_index_`: `int`
  * `value_`: `Eigen::VectorXd`

### **NeumannBC**
* **Description**: Applies a known flux to a node, modifying the force vector F.
* **Public Functions**:
  * `NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value)`
  * `apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override`
* **Private Members**:
  * `equation_index_`: `int`
  * `flux_value_`: `Eigen::VectorXd`

### **CauchyBC**
* **Description**: Models mixed/convective boundaries, modifying both K and F.
* **Public Functions**:
  * `CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf)`
  * `apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override`
* **Private Members**:
  * `equation_index_`: `int`
  * `h_`: `Eigen::VectorXd`
  * `T_inf_`: `Eigen::VectorXd`

### **DOFManager**
* **Description**: The Degree of Freedom Manager. It maps the physics variables at each node to a unique index in the global system of equations.
* **Public Functions**:
  * `DOFManager(Mesh& mesh)`
  * `registerVariable(const std::string& var_name)`
  * `build()`
  * `getEquationIndex(int node_id, const std::string& var_name) const`
  * `getNumEquations() const`: returns `size_t`
  * `getVariableNames() const`: returns `const std::vector<std::string>&`
* **Private Members**:
  * `mesh_`: `Mesh&`
  * `variable_names_`: `std::vector<std::string>`
  * `dof_map_`: `std::map<std::pair<int, int>, int>`
  * `num_equations_`: `size_t`

### **Element**
* **Description**: An abstract base class for all mesh elements.
* **Public Functions**:
  * `Element(int id)`
  * `getId() const`: returns `int`
  * `getNodes() const`: returns `const std::vector<Node*>&`
  * `addNode(Node* node)`
  * `getNumNodes() const`: returns `size_t` (pure virtual)
  * `getTypeName() const`: returns `const char*` (pure virtual)
* **Protected Members**:
  * `id_`: `int`
  * `nodes_`: `std::vector<Node*>`
* **Derived Classes**:
  * `LineElement`
  * `TriElement`
  * `TetElement`

### **LineElement**
* **Description**: A concrete implementation for a 1D line element with 2 nodes.
* **Public Functions**:
  * `LineElement(int id)`
  * `getNumNodes() const override`: returns `size_t`
  * `getTypeName() const override`: returns `const char*`
  * `getLength() const`: returns `double`

### **TriElement**
* **Description**: Represents a 2D, 3-node, linear triangular element.
* **Public Functions**:
  * `TriElement(int id)`
  * `getNumNodes() const override`: returns `size_t`
  * `getTypeName() const override`: returns `const char*`
  * `getArea() const`: returns `double`
  * `getBMatrix() const`: returns `Eigen::Matrix<double, 2, 3>`

### **TetElement**
* **Description**: Represents a 3D, 4-node, linear tetrahedral element.
* **Public Functions**:
  * `TetElement(int id)`
  * `getNumNodes() const override`: returns `size_t`
  * `getTypeName() const override`: returns `const char*`
  * `getVolume() const`: returns `double`
  * `getBMatrix() const`: returns `Eigen::Matrix<double, 3, 4>`

### **LinearSolver**
* **Description**: A simple, static wrapper class around the Eigen library's sparse solver (`SparseLU`). It provides a single function `solve(A, b, x)`.
* **Public Functions**:
  * `solve(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& b, Eigen::MatrixXd& x)`: Solves the system Ax = b and returns x.

### **Material**
* **Description**: Encapsulates the physical properties of a material. This class is designed to handle both constant and temperature-dependent properties.
* **Public Functions**:
  * `Material(const std::string& name)`
  * `setProperty(const std::string& prop_name, double value)`
  * `setTempDependentProperty(const std::string& prop_name, const std::map<std::string, double>& params)`
  * `getProperty(const std::string& prop_name) const`: returns `double`
  * `getProperty(const std::string& prop_name, double temperature) const`: returns `double`
  * `getName() const`: returns `const std::string&`
* **Private Members**:
  * `name_`: `std::string`
  * `properties_`: `std::map<std::string, std::any>`

### **Mesh**
* **Description**: A container for all nodes and elements in the simulation.
* **Public Functions**:
  * `~Mesh()`
  * `addNode(Node* node)`
  * `addElement(Element* element)`
  * `getNode(int id) const`: returns `Node*`
  * `getElement(int id) const`: returns `Element*`
  * `getElements() const`: returns `const std::vector<Element*>&`
  * `getNodes() const`: returns `const std::vector<Node*>&`
  * `create_uniform_1d_mesh(double length, int num_elements)`: returns `Mesh*`
  * `create_uniform_2d_mesh(double width, double height, int nx, int ny)`: returns `Mesh*`
  * `create_uniform_3d_mesh(double width, double height, double depth, int nx, int ny, int nz)`: returns `Mesh*`
* **Private Members**:
  * `nodes_`: `std::vector<Node*>`
  * `elements_`: `std::vector<Element*>`
  * `node_map_`: `std::map<int, Node*>`
  * `element_map_`: `std::map<int, Element*>`

### **Node**
* **Description**: Represents a point in space with an ID and coordinates.
* **Public Functions**:
  * `Node(int id, double x, double y = 0.0, double z = 0.0)`
  * `getId() const`: returns `int`
  * `getCoords() const`: returns `const std::vector<double>&`
* **Private Members**:
  * `id_`: `int`
  * `coords_`: `std::vector<double>`

### **Problem**
* **Description**: The main orchestrator class. It owns all the components of a simulation (mesh, fields, BCs) and controls the high-level workflow.
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
* **Private Members**:
  * `mesh_`: `std::unique_ptr<Mesh>`
  * `dof_manager_`: `std::unique_ptr<DOFManager>`
  * `fields_`: `std::vector<std::unique_ptr<Physics::PhysicsField>>`
  * `max_iterations_`: `int`
  * `convergence_tolerance_`: `double`
  * `time_step_`: `double`
  * `total_time_`: `double`