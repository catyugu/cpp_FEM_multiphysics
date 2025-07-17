# **PhysicsField Namespace**

This namespace contains the specific implementations for different physical simulations. Its core design principle is to abstract the mathematical formulation of a physical phenomenon away from the main solver and mesh components.

The most powerful feature of this namespace is its implementation of a **P-Refinement** (or P-enrichment) strategy, allowing for higher-order mathematical approximations on simple, geometrically linear meshes.

---
## **Classes**

### **PhysicsField**
* **Description**: This is an **abstract base class** for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module. It provides the framework for assembling system matrices and applying boundary conditions and sources.
* **Key Concept: P-Refinement Strategy**:
  * The physical `Mesh` is always composed of simple, geometrically linear elements (e.g., a triangle in the mesh data structure always has exactly 3 vertex nodes).
  * The accuracy of the simulation can be improved by increasing the **mathematical order** of the shape functions used to approximate the solution over these simple elements. This is controlled by the `element_order_` property.
  * Setting `element_order_` to `2` instructs the `assemble` method to behave as if the element were quadratic. It achieves this by creating "virtual" nodes on the midpoints of the element's edges and using 2nd-order shape functions from `Utils::ShapeFunctions`. This entire process is handled internally within the `assemble` method and does not require a different mesh file.
* **Public Functions**:
  * `~PhysicsField()`: Virtual destructor.
  * `getName() const`: Returns the name of the physics field (e.g., "Heat Transfer 2D").
  * `getVariableName() const`: Returns the name of the primary variable (e.g., "Temperature").
  * `getMaterial() const`: Returns a constant reference to the `Core::Material` object.
  * `getDimension() const`: Returns the physical dimension of the field (1, 2, or 3). This is a **new virtual function**.
  * `setElementOrder(int order)` and `getElementOrder() const`: Sets and gets the mathematical order of the approximation (e.g., 1 for linear, 2 for quadratic).
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager)`: Sets up the field's internal matrices based on the total number of equations determined by the `DOFManager`.
  * `assemble()`: Assembles the global stiffness (`K`) and mass (`M`) matrices. This is the core numerical method. It now uses the `FEValues` object created by the `Element` to simplify calculations significantly. It iterates through each mesh element and performs the following steps:
    1.  Calls `elem->create_fe_values(element_order_)` to get the pre-calculated `FEValues` object for the element.
    2.  Gathers the Degrees of Freedom (DOFs) for both the real vertex nodes and the "virtual" higher-order nodes (e.g., edge midpoints) using the `DOFManager` and the `get_element_dofs` helper.
    3.  Loops over quadrature points using `fe_values->reinit(q_p)` and retrieves pre-calculated `N`, `B` (gradients in real coordinates), and `detJ_x_w` (Jacobian determinant times quadrature weight).
    4.  Computes local stiffness and/or mass matrices using these values.
    5.  Adds the local matrices' contributions to the correct locations in the global system matrices.
  * `addBC(std::unique_ptr<Core::BoundaryCondition> bc)`: Adds a boundary condition to a list.
  * `applyBCs()`: Applies all stored boundary conditions to the system. This method intelligently uses a consolidation map to ensure each degree of freedom is constrained by only one `DirichletBC` object, preventing redundant constraints.
  * `addSource(std::unique_ptr<Core::SourceTerm> source)` and `applySources()`: Manages the application of domain source terms to the right-hand-side vector `F`. `applySources()` now clears `F_` before applying sources to prevent accumulation.
  * `setInitialConditions(...)`: Sets initial values for the field. Overloaded for constant values and functions.
  * `getStiffnessMatrix()`, `getMassMatrix()`, `getRHS()`, `getSolution()`: Accessors for matrices and vectors.
  * `getPreviousSolution() const`: Returns the solution from the previous time step (for transient problems).
  * `updatePreviousSolution()`: Updates `U_prev_` with the current solution `U_`.
  * `getMesh() const`: Returns a pointer to the mesh.
  * `getDofManager() const`: Returns a pointer to the DOFManager.
  * `removeBCsByTag(const std::string& tag)`: Removes boundary conditions associated with a specific tag.
  * `removeSourcesByTag(const std::string& tag)`: Removes source terms associated with a specific tag.
  * `enable()`, `disable()`, `isEnabled()`: Control whether a field is active in the simulation.
* **Protected Members**:
  * `mesh_`, `dof_manager_`, `element_order_`: Pointers to core components and the mathematical order.
  * `K_`, `M_`, `F_`, `U_`, `U_prev_`: System matrices and vectors.
  * `bcs_`, `source_terms_`: Vectors holding the boundary conditions and source terms.
  * `get_element_dofs(Core::Element* elem) const`: A helper function to collect the global DOF indices for an element, considering both vertex and higher-order (edge/face/volume) nodes based on the `element_order_`. This function now handles the canonical ordering of higher-order nodes (midpoints for order 2) for Line, Tri, and Tet elements.

---
### **Derived Physics Classes**
* **`Current1D/2D/3D`**
* **`Heat1D/2D/3D`**
* **`Magnetic1D/2D/3D`**

* **Description**: These are the concrete implementations of `PhysicsField` for different physical phenomena. They each implement the `assemble` method with the specific material properties (e.g., electrical conductivity, thermal conductivity, magnetic permeability) and the correct PDE formulation for their respective physics. All of these classes now correctly implement the advanced assembly logic required to support higher-order element approximations by leveraging the `FEValues` object.