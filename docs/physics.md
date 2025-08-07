# **Physics Namespace**

This namespace contains the specific implementations for different physical simulations. Its core design principle is to abstract the mathematical formulation of a physical phenomenon away from the main solver and mesh components.

A key feature of this namespace is its implementation of a **P-Refinement** strategy, allowing for higher-order mathematical approximations on simple, geometrically linear meshes.

---
## **Classes**

### **PhysicsField**
* **Description**: This is an **abstract base class** for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module.
* **Key Features**:
  * **P-Refinement**: Accuracy is controlled by the `element_order_` property. Setting this to `2` or higher instructs the assembly process to use higher-order shape functions over the simple geometric elements.
  * **Vector Field Support**: The `getNumComponents()` virtual function allows fields to declare themselves as scalar (`return 1;`) or vector (`return 3`), a feature used by `Magnetic3D`.
  * **Coupling Support**: Contains a separate right-hand-side vector `F_coupling_` to cleanly manage source terms arising from multiphysics interactions.
  * **Material Handling**: Physics fields no longer store a single material instance. Instead, they retrieve material properties on a per-element basis through the `Problem` class, allowing for simulations with multiple materials.
* **Public Functions**:
  * `getName() const`, `getVariableName() const`, `getDimension() const`: Pure virtual functions to define the field's identity.
  * `getMaterial(const Core::Element* elem) const`: Pure virtual function to retrieve material properties for a specific element. This enables support for spatially varying material properties within a single simulation.
  * `getNumComponents() const`: Virtual function to specify if the field is scalar (1 component) or vector (N components). Defaults to 1.
  * `setElementOrder(int order)`: Sets the mathematical order of the approximation (e.g., 1 for linear, 2 for quadratic).
  * `setup(...)`: Initializes the field's matrices and vectors based on the DOF map.
  * `assemble(const PhysicsField *coupled_field = nullptr)`: Assembles the global stiffness (`K`) and mass (`M`) matrices. This is the core numerical method. It iterates through elements and uses the `ReferenceElementCache` and `FEValues` calculator to efficiently compute local matrices. Material properties are now retrieved per element rather than using a single global material.
  * `addBC(...)`, `applyBCs()`: Manages boundary conditions. `applyBCs()` intelligently consolidates all Dirichlet BCs to prevent redundant constraints on shared nodes.
  * `addSource(...)`, `applySources()`: Manages domain source terms.
  * Accessors for matrices (`K_`, `M_`), vectors (`F_`, `U_`, `U_prev_`, `F_coupling_`), mesh, and DOF manager.
  * `getElementDofs(Core::Element* elem) const`: A helper function to collect all global DOF indices for an element, correctly handling both vertex and higher-order (edge) nodes based on `element_order_` and `num_components_`.

---### **Derived Physics Classes**
* **`Current1D/2D/3D`**, **`Heat1D/2D/3D`**, **`Magnetic1D/2D/3D`**
  * These are concrete implementations for scalar fields (Voltage, Temperature, Magnetic Potential). They override `getNumComponents()` to return `1`. Their `assemble` methods implement the weak form of the corresponding scalar PDE (e.g., Laplace or Poisson equation). Material properties are now retrieved per element rather than using a constructor-provided single material instance.
  * **Material Handling**: These classes no longer accept a material in their constructors. Instead, they retrieve material properties from the `Problem` class on a per-element basis using the `getMaterial(const Core::Element* elem)` method during assembly.

* **`Magnetic3D`** (Magnetostatics)
  * **Description**: A concrete implementation for solving 3D magnetostatic problems based on the magnetic vector potential **A**.
  * **Variable**: `MagneticVectorPotential`
  * **Components**: Overrides `getNumComponents()` to return **3**.
  * **`assemble()` Method**: This method is the core of the magnetostatics solver. It builds the local stiffness matrix for each element by implementing the weak form of the magnetostatic equation, $\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) = \mathbf{J}$. This involves using the **curl of the vector shape functions**, which is constructed from the shape function gradients (`∇N`) provided by `FEValues`. Material properties (magnetic permeability) are retrieved per element.