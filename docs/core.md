# **Core Namespace**

This namespace contains the fundamental, high-level components that orchestrate the simulation framework. These classes manage the overall problem setup, degrees of freedom, material properties, boundary conditions, and domain sources.

---
## **Classes**

### **Problem**

* **Description**: The main orchestrator class. It owns all the components of a simulation (mesh, fields, solver) and controls the high-level workflow.
* **Public Functions**:
  * `Problem(std::unique_ptr<Mesh> mesh)`
  * `~Problem()`
  * `addField(std::unique_ptr<Physics::PhysicsField> field)`: Adds a physics field to the problem. It correctly reads the number of components from the field (scalar or vector) and passes this information to the `DOFManager`.
  * `addPostProcessor(std::unique_ptr<Post::PostProcessor> post_processor)`: Adds a post-processing calculator to the problem (e.g., for calculating heat flux).
  * `setup()`: Collects the element order for each physics field and passes it to the `DOFManager::build()` method, ensuring correct DOF mapping for higher-order elements. It also sets up the appropriate solver via a factory.
  * `setIterativeSolverParameters(int max_iter, double tol)`: Sets parameters for iterative solvers.
  * `setTimeSteppingParameters(double time_step, double total_time)`: Sets parameters for transient simulations.
  * `setLinearSolverType(Solver::SolverType type)`: Specifies the linear solver to be used (LU or BiCGSTAB).
  * `solveSteadyState()` and `solveTransient()`: Executes the appropriate solver strategy and subsequently runs any registered post-processors.
  * `exportResults(const std::string& filename) const`
  * Accessors like `getMesh()`, `getDofManager()`, `getField()`, `getPostProcessingResults()`, etc.
* **Private Members**:
  * `mesh_`, `dof_manager_`, `solver_`, `fields_`: Pointers to the core components.
  * `post_processors_`, `post_processing_results_`: Manages post-processing tasks.
  * `coupling_manager_`: Manages interactions between different physics fields.

### **DOFManager**

* **Description**: The Degree of Freedom (DOF) Manager. It intelligently maps physics variables to unique indices in the global system of equations. It is designed to handle both **scalar and vector variables**, as well as higher-order element formulations.
* **Public Functions**:
  * `DOFManager(Mesh& mesh)`
  * `registerVariable(const std::string& var_name, int num_components = 1)`: Registers a variable with a specified number of components. For a scalar field, `num_components` is 1. For a 3D vector field, it's 3.
  * `build(const std::map<std::string, int>& field_orders)`: Builds the complete DOF map. It correctly allocates a contiguous block of equation indices for each component of a vector variable at every node (both vertex and higher-order).
  * `getEquationIndex(int node_id, const std::string& var_name) const`: Gets the global equation index for a variable at a vertex. For a vector field, this is the index of the *first component* (e.g., Ax).
  * `getEdgeEquationIndex(std::vector<int> node_ids, const std::string& var_name) const`: Gets the global equation index for a higher-order edge DOF.

### **Material**

* **Description**: Encapsulates the physical properties of a material. This class is designed to handle both constant properties and properties that are dependent on other field values (e.g., temperature-dependent electrical conductivity) by using `std::function`.

### **BoundaryCondition**

* **Description**: An abstract base class for applying conditions on the boundary of the domain.
* **Derived Classes**:
  * `DirichletBC`: Has two constructors: one for standard vertex nodes (taking a `node_id` and `var_name`) and a second, more direct constructor for any DOF (taking an `equation_index` directly). This second constructor is crucial for applying constraints to individual components of a vector field or to higher-order "virtual" nodes. It also provides a static factory `create` method for conveniently generating BCs based on a geometric predicate.
  * `NeumannBC`: Applies a flux condition to a specific degree of freedom.
  * `CauchyBC`: Applies a mixed boundary condition (e.g., for convection) to a specific degree of freedom.

### **SourceTerm**

* **Description**: An abstract base class for applying sources within the domain (e.g., a heat source or current density).
* **Key Update**: The `apply` method signature now includes `int element_order` to correctly handle the assembly for higher-order elements.
* **Derived Classes**:
  * `VolumetricSource`: Applies a total power source over the volume of an element.
  * `PrescribedCurrentDensity`: Applies a vector-valued electric current density **J** to an element, used as the source term for the `Magnetic3D` field.