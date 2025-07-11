# **PhysicsField Namespace**

This namespace contains the specific implementations for different physical simulations.

---
## **Classes**

### **PhysicsField**
* **Description**: This is an abstract base class for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module.
* **Public Functions**:
    * `~PhysicsField()`
    * `getName() const`: returns `const char*` (pure virtual)
    * `getVariableName() const`: returns `const char*` (pure virtual)
    * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager)` (pure virtual)
    * `assemble()` (pure virtual)
    * `addBC(std::unique_ptr<Core::BoundaryCondition> bc)`
    * `applyBCs()`
    * `getBCs() const`: returns `const std::vector<std::unique_ptr<Core::BoundaryCondition>>&`
    * `updatePreviousSolution()`
    * `getPreviousSolution() const`: returns `const Eigen::VectorXd&`
    * `setInitialConditions(double initial_value)`
    * `setInitialConditions(std::function<F> initial_conditions)`
    * `getStiffnessMatrix()`: returns `Eigen::SparseMatrix<double>&`
    * `getMassMatrix()`: returns `Eigen::SparseMatrix<double>&`
    * `getRHSVector()`: returns `Eigen::VectorXd&`
    * `getSolution()`: returns `Eigen::VectorXd&`
    * `getStiffnessMatrix() const`: returns `const Eigen::SparseMatrix<double>&`
    * `getMassMatrix() const`: returns `const Eigen::SparseMatrix<double>&`
    * `getRHSVector() const`: returns `const Eigen::VectorXd&`
    * `getSolution() const`: returns `const Eigen::VectorXd&`
* **Protected Members**:
    * `mesh_`: `Core::Mesh*`
    * `dof_manager_`: `Core::DOFManager*`
    * `K_`: `Eigen::SparseMatrix<double>`
    * `M_`: `Eigen::SparseMatrix<double>`
    * `F_`: `Eigen::VectorXd`
    * `U_`: `Eigen::VectorXd`
    * `U_prev_`: `Eigen::VectorXd`
    * `bcs_`: `std::vector<std::unique_ptr<Core::BoundaryCondition>>`

### **Current1D**
* **Description**: Implements the `assemble` method for 1D Electromagnetics.
* **Public Functions**:
    * `Current1D(const Core::Material& material)`
    * `getName() const override`: returns `const char*`
    * `getVariableName() const override`: returns `const char*`
    * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`
    * `assemble() override`
    * `calculateJouleHeat() const`: returns `std::vector<double>`
    * `setCoupledHeatField(const PhysicsField* heat_field)`
* **Private Members**:
    * `material_`: `const Core::Material&`
    * `heat_field_`: `const PhysicsField*`

### **Current2D**
* **Description**: Implements the `assemble` method for 2D Electromagnetics.
* **Public Functions**:
    * `Current2D(const Core::Material& material)`
    * `getName() const override`: returns `const char*`
    * `getVariableName() const override`: returns `const char*`
    * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`
    * `assemble() override`
    * `calculateJouleHeat() const`: returns `std::vector<double>`
    * `setCoupledHeatField(const PhysicsField* heat_field)`
* **Private Members**:
    * `material_`: `const Core::Material&`
    * `heat_field_`: `const PhysicsField*`

### **Heat1D**
* **Description**: Implements the `assemble` method for 1D Heat Transfer.
* **Public Functions**:
    * `Heat1D(const Core::Material& material)`
    * `getName() const override`: returns `const char*`
    * `getVariableName() const override`: returns `const char*`
    * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`
    * `assemble() override`
    * `setVolumetricHeatSource(const std::vector<double>& source)`
* **Private Members**:
    * `material_`: `const Core::Material&`
    * `k_`: `double`
    * `rho_`: `double`
    * `cp_`: `double`
    * `volumetric_heat_source_`: `std::vector<double>`

### **Heat2D**
* **Description**: Implements the `assemble` method for 2D Heat Transfer.
* **Public Functions**:
    * `Heat2D(const Core::Material& material)`
    * `getName() const override`: returns `const char*`
    * `getVariableName() const override`: returns `const char*`
    * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`
    * `assemble() override`
    * `setVolumetricHeatSource(const std::vector<double>& source)`
* **Private Members**:
    * `material_`: `const Core::Material&`
    * `k_`: `double`
    * `rho_`: `double`
    * `cp_`: `double`
    * `volumetric_heat_source_`: `std::vector<double>`

### **Heat3D**
* **Description**: Implements the `assemble` method for 3D Heat Transfer.
* **Public Functions**:
    * `Heat3D(const Core::Material& material)`
    * `getName() const override`: returns `const char*`
    * `getVariableName() const override`: returns `const char*`
    * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override`
    * `assemble() override`
* **Private Members**:
    * `material_`: `const Core::Material&`
    * `k_`: `double`

