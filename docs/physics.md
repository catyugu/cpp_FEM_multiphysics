# **PhysicsField Namespace**

This namespace contains the specific implementations for different physical simulations.

---
## **Classes**

### **PhysicsField**
* **Description**: This is an **abstract base class** for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module.
* **Public Functions**:
  * `~PhysicsField()`: Virtual destructor.
  * `getName() const`: Returns the name of the physics field (e.g., "Heat Transfer 2D") as a `const char*`. This is a pure virtual function.
  * `getVariableName() const`: Returns the name of the primary variable for this physics field (e.g., "Temperature") as a `const char*`. This is a pure virtual function.
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager)`: Sets up the physics field with the given mesh and DOF manager. This is a pure virtual function.
  * `assemble()`: Assembles the stiffness matrix (K), mass matrix (M), and force vector (F) for the physics field. This is a pure virtual function.
  * `addBC(std::unique_ptr<Core::BoundaryCondition> bc)`: Adds a boundary condition to the physics field.
  * `applyBCs()`: Applies all added boundary conditions to the system matrices.
  * `getBCs() const`: Returns a constant reference to the vector of boundary conditions (`const std::vector<std::unique_ptr<Core::BoundaryCondition>>&`).
  * `updatePreviousSolution()`: Copies the current solution (`U_`) to the previous solution (`U_prev_`).
  * `getPreviousSolution() const`: Returns a constant reference to the solution from the previous time step (`const Eigen::MatrixXd&`).
  * `setInitialConditions(double initial_value)`: Sets a uniform initial value for the entire field.
  * `setInitialConditions(std::function<F> initial_conditions)`: Sets initial conditions using a function that takes a node pointer and returns a value.
  * `getStiffnessMatrix()`: Returns a reference to the stiffness matrix (`Eigen::SparseMatrix<double>&`).
  * `getMassMatrix()`: Returns a reference to the mass matrix (`Eigen::SparseMatrix<double>&`).
  * `getRHS()`: Returns a reference to the right-hand side (force) vector (`Eigen::MatrixXd&`).
  * `getSolution()`: Returns a reference to the solution vector (`Eigen::MatrixXd&`).
  * `getStiffnessMatrix() const`: Returns a constant reference to the stiffness matrix (`const Eigen::SparseMatrix<double>&`).
  * `getMassMatrix() const`: Returns a constant reference to the mass matrix (`const Eigen::SparseMatrix<double>&`).
  * `getRHS() const`: Returns a constant reference to the right-hand side vector (`const Eigen::MatrixXd&`).
  * `getSolution() const`: Returns a constant reference to the solution vector (`const Eigen::MatrixXd&`).
  * `enable()`: Enables the physics field.
  * `disable()`: Disables the physics field.
  * `isEnabled() const`: Returns true if the physics field is enabled.
* **Protected Members**:
  * `mesh_`: A pointer to the `Core::Mesh` object (`Core::Mesh*`).
  * `dof_manager_`: A pointer to the `Core::DOFManager` object (`Core::DOFManager*`).
  * `enabled`: A boolean indicating if the field is active (`bool`).
  * `K_`: The stiffness matrix (`Eigen::SparseMatrix<double>`).
  * `M_`: The mass matrix (`Eigen::SparseMatrix<double>`).
  * `F_`: The right-hand side (force) vector (`Eigen::MatrixXd`).
  * `U_`: The current solution vector (`Eigen::MatrixXd`).
  * `U_prev_`: The solution vector from the previous time step (`Eigen::MatrixXd`).
  * `bcs_`: A vector of unique pointers to the boundary conditions (`std::vector<std::unique_ptr<Core::BoundaryCondition>>`).

### **Current1D**
* **Description**: Implements the `assemble` method for 1D Electromagnetics.
* **Public Functions**:
  * `Current1D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Electromagnetics 1D".
  * `getVariableName() const override`: Returns "Voltage".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 1D current field.
  * `assemble() override`: Assembles the system matrices for the 1D current field.
  * `calculateJouleHeat() const`: Calculates the Joule heating for each element and returns it as a `std::vector<double>`.
  * `setCoupledHeatField(const PhysicsField* heat_field)`: Links this EMag field to a heat field for temperature-dependent calculations.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).
  * `heat_field_`: A pointer to the coupled `PhysicsField` for heat (`const PhysicsField*`).

### **Current2D**
* **Description**: Implements the `assemble` method for 2D Electromagnetics.
* **Public Functions**:
  * `Current2D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Electromagnetics 2D".
  * `getVariableName() const override`: Returns "Voltage".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 2D current field.
  * `assemble() override`: Assembles the system matrices for the 2D current field.
  * `calculateJouleHeat() const`: Calculates the Joule heating for each element and returns it as a `std::vector<double>`.
  * `setCoupledHeatField(const PhysicsField* heat_field)`: Links this EMag field to a heat field for temperature-dependent calculations.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).
  * `heat_field_`: A pointer to the coupled `PhysicsField` for heat (`const PhysicsField*`).

### **Heat1D**
* **Description**: Implements the `assemble` method for 1D Heat Transfer.
* **Public Functions**:
  * `Heat1D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Heat Transfer 1D".
  * `getVariableName() const override`: Returns "Temperature".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 1D heat field.
  * `assemble() override`: Assembles the system matrices for the 1D heat field.
  * `setVolumetricHeatSource(const std::vector<double>& source)`: Sets the volumetric heat source for each element.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).
  * `k_`: Thermal conductivity (`double`).
  * `rho_`: Density (`double`).
  * `cp_`: Specific heat capacity (`double`).
  * `volumetric_heat_source_`: A vector containing the volumetric heat source for each element (`std::vector<double>`).

### **Heat2D**
* **Description**: Implements the `assemble` method for 2D Heat Transfer.
* **Public Functions**:
  * `Heat2D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Heat Transfer 2D".
  * `getVariableName() const override`: Returns "Temperature".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 2D heat field.
  * `assemble() override`: Assembles the system matrices for the 2D heat field.
  * `setVolumetricHeatSource(const std::vector<double>& source)`: Sets the volumetric heat source for each element.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).
  * `k_`: Thermal conductivity (`double`).
  * `rho_`: Density (`double`).
  * `cp_`: Specific heat capacity (`double`).
  * `volumetric_heat_source_`: A vector containing the volumetric heat source for each element (`std::vector<double>`).

### **Heat3D**
* **Description**: Implements the `assemble` method for 3D Heat Transfer.
* **Public Functions**:
  * `Heat3D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Heat Transfer 3D".
  * `getVariableName() const override`: Returns "Temperature".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 3D heat field.
  * `assemble() override`: Assembles the system matrices for the 3D heat field.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).
  * `k_`: Thermal conductivity (`double`).

Of course. Here is the updated documentation reflecting the addition of the `Magnetic1D` physics module.

# **PhysicsField Namespace**

This namespace contains the specific implementations for different physical simulations.

---
## **Classes**
... (existing classes: PhysicsField, Current1D, Current2D, Heat1D, Heat2D, Heat3D) ...

### **Magnetic1D**
* **Description**: Implements the `assemble` method for 1D Magnetostatics using the magnetic vector potential, `A_z`, as the primary unknown.
* **Public Functions**:
  * `Magnetic1D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Magnetic Field 1D".
  * `getVariableName() const override`: Returns "MagneticPotential".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 1D magnetic field.
  * `assemble() override`: Assembles the system matrices for the 1D magnetic field.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).

### **Magnetic2D**
* **Description**: Implements the `assemble` method for 2D Magnetostatics. The formulation uses the magnetic vector potential, `A_z`, as the primary unknown.
* **Public Functions**:
  * `Magnetic2D(const Core::Material& material)`: Constructor that takes a material object.
  * `getName() const override`: Returns "Magnetic Field 2D".
  * `getVariableName() const override`: Returns "MagneticPotential".
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`: Sets up the 2D magnetic field.
  * `assemble() override`: Assembles the system matrices for the 2D magnetic field based on the equation ∇ ⋅ ( (1/μ) ∇A_z ) = -J_z.
* **Private Members**:
  * `material_`: A constant reference to the `Core::Material` object (`const Core::Material&`).
