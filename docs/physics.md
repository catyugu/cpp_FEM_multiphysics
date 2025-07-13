# **PhysicsField Namespace**

This namespace contains the specific implementations for different physical simulations.

---
## **Classes**

### **PhysicsField**
* **Description**: This is an **abstract base class** for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module.
* **Public Functions**:
  * `~PhysicsField()`: Virtual destructor.
  * `getName() const`: Returns the name of the physics field (e.g., "Heat Transfer 2D"). Pure virtual.
  * `getVariableName() const`: Returns the name of the primary variable (e.g., "Temperature"). Pure virtual.
  * `getMaterial() const`: Returns a constant reference to the `Core::Material` object. Pure virtual.
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager)`: Sets up the field. Pure virtual.
  * `assemble()`: Assembles the system matrices (K, M). Pure virtual.
  * `addBC(std::unique_ptr<Core::BoundaryCondition> bc)`: Adds a boundary condition.
  * `removeBCsByTag(const std::string& tag)`: Removes all BCs with a matching tag.
  * `applyBCs()`: Applies all added boundary conditions to the system matrices.
  * `addSource(std::unique_ptr<Core::SourceTerm> source)`: Adds a domain source term.
  * `removeSourcesByTag(const std::string& tag)`: Removes all sources with a matching tag.
  * `applySources()`: Assembles the force vector `F` from all added source terms.
  * `setInitialConditions(...)`: Sets initial values for the field.
  * `updatePreviousSolution()`: Copies the current solution to the previous time step.
  * `getStiffnessMatrix()`, `getMassMatrix()`, `getRHS()`, `getSolution()`: Accessors for matrices and vectors.
* **Protected Members**:
  * `mesh_`, `dof_manager_`, `enabled`, `element_order_`.
  * `K_`, `M_`, `F_`, `U_`, `U_prev_`: System matrices and vectors.
  * `bcs_`: Vector of boundary conditions.
  * `source_terms_`: Vector of source terms.

### **Current1D / Current2D / Current3D**
* **Description**: Implements the `assemble` method for 1D and 2D Electromagnetics ("Voltage"). The classes are now purely responsible for assembling their own systems and are no longer aware of thermal coupling.

### **Heat1D / Heat2D / Heat3D**
* **Description**: Implements the `assemble` method for 1D and 2D Heat Transfer ("Temperature"). The classes no longer contain specialized methods for volumetric heat sources; all sources are now handled by the generic `SourceTerm` system.

### **Magnetic1D / Magnetic2D / Magnetic3D**
* **Description**: Implements the `assemble` method for 1D and 2D Magnetostatics ("MagneticPotential").