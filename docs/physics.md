Of course. Here is the complete and detailed documentation for the `physics.md` file, updated to reflect the final, robust architecture for handling higher-order elements.

***

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
  * `getDimension() const`: Returns the physical dimension of the field (1, 2, or 3).
  * `setElementOrder(int order)` and `getElementOrder() const`: Sets and gets the mathematical order of the approximation (e.g., 1 for linear, 2 for quadratic).
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager)`: Sets up the field's internal matrices based on the total number of equations determined by the `DOFManager`.
  * `assemble()`: Assembles the global stiffness (`K`) and mass (`M`) matrices. This is the core numerical method. It iterates through each mesh element and performs the following steps:
    1.  Determines the number of nodes for the "virtual" higher-order element based on `element_order_`.
    2.  Gathers the Degrees of Freedom (DOFs) for both the real vertex nodes and the "virtual" higher-order nodes (e.g., edge midpoints) using the `DOFManager`.
    3.  Calculates the local element stiffness and mass matrices using numerical integration (`Utils::Quadrature`) and the appropriate higher-order shape functions (`Utils::ShapeFunctions`).
    4.  Adds the local matrices' contributions to the correct locations in the global system matrices.
  * `addBC(std::unique_ptr<Core::BoundaryCondition> bc)`: Adds a boundary condition to a list.
  * `applyBCs()`: Applies all stored boundary conditions to the system. This method is now intelligent, using a consolidation map to ensure that each degree of freedom is constrained by only one `DirichletBC` object, even if multiple are defined for the same node in the test or main setup. This prevents solver failures from redundant constraints.
  * `addSource(std::unique_ptr<Core::SourceTerm> source)` and `applySources()`: Manages the application of domain source terms to the right-hand-side vector `F`.
  * `setInitialConditions(...)`: Sets initial values for the field.
  * `getStiffnessMatrix()`, `getMassMatrix()`, `getRHS()`, `getSolution()`: Accessors for matrices and vectors.

* **Protected Members**:
  * `mesh_`, `dof_manager_`, `element_order_`: Pointers to core components and the mathematical order.
  * `K_`, `M_`, `F_`, `U_`, `U_prev_`: System matrices and vectors.
  * `bcs_`, `source_terms_`: Vectors holding the boundary conditions and source terms.

---
### **Derived Physics Classes**
* **`Current1D/2D/3D`**
* **`Heat1D/2D/3D`**
* **`Magnetic1D/2D/3D`**

* **Description**: These are the concrete implementations of `PhysicsField` for different physical phenomena. They each implement the `assemble` method with the specific material properties (e.g., electrical conductivity, thermal conductivity, magnetic permeability) and the correct PDE formulation for their respective physics. All of these classes now correctly implement the advanced assembly logic required to support higher-order element approximations.