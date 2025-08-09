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
  * `assemble()`: Assembles the global stiffness (`K`) and mass (`M`) matrices. This is the core numerical method. The process is:
    1. Iterate through each element in the mesh.
    2. For each element, iterate through its quadrature points.
    3. At each quadrature point:
        a. Interpolate the required field variables (e.g., Temperature) to that specific point using `Utils::InterpolationUtilities`.
        b. Retrieve the material for the element and evaluate its properties (e.g., conductivity) at that point using the interpolated variables via `material.getPropertyAtQuadraturePoint()`.
        c. Calculate the contribution to the element's local stiffness/mass matrix using the local material properties and shape function data from `FEValues`.
    4. Add the completed local matrix to the global system matrix.
    This per-quadrature-point approach ensures high accuracy for non-linear and coupled problems.
  * `addBC(...)`, `applyBCs()`: Manages boundary conditions. `applyBCs()` intelligently consolidates all Dirichlet BCs to prevent redundant constraints on shared nodes.
  * `addSource(...)`, `applySources()`: Manages domain source terms.
  * Accessors for matrices (`K_`, `M_`), vectors (`F_`, `U_`, `U_prev_`, `F_coupling_`), mesh, and DOF manager.
  * `getElementDofs(Core::Element* elem) const`: A helper function to collect all global DOF indices for an element, correctly handling both vertex and higher-order (edge) nodes based on `element_order_` and `num_components_`.

---### **Derived Physics Classes**
* **`Current1D/2D/3D`**, **`Heat1D/2D/3D`**, **`Magnetic1D/2D/3D`**
  * These are concrete implementations for scalar fields (Voltage, Temperature, Magnetic Potential).
  * **Material Handling**: Their `assemble` methods now follow the per-quadrature-point evaluation workflow described above to accurately capture non-linear material behavior. They no longer use a single property value for the entire element.

* **`Magnetic3D`** (Magnetostatics)
  * **Description**: A concrete implementation for solving 3D magnetostatic problems based on the magnetic vector potential **A**.
  * **Variable**: `MagneticVectorPotential`
  * **Components**: Overrides `getNumComponents()` to return **3**.
  * **`assemble()` Method**: This method is the core of the magnetostatics solver. It builds the local stiffness matrix for each element by implementing the weak form of the magnetostatic equation, $\nabla \times \left( \frac{1}{\mu} \nabla \times \mathbf{A} \right) = \mathbf{J}$. It uses the curl of the vector shape functions and evaluates the magnetic permeability (`mu`) at each quadrature point to handle non-linear magnetic materials.
