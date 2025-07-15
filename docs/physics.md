# **PhysicsField Namespace**

This namespace contains the specific implementations for different physical simulations.

---
## **Classes**

### **PhysicsField**
* **Description**: This is an **abstract base class** for all physics implementations. It defines the common interface that the `Problem` class uses to interact with any physics module. It now supports **P-Refinement** (higher-order elements).
* **Key Concept: P-Refinement**:
  * The mesh is always composed of simple, geometrically linear elements (e.g., 3-node triangles).
  * The accuracy of the simulation can be improved by increasing the mathematical order of the shape functions used to approximate the solution over these simple elements.
  * Setting `element_order_` to `2` instructs the `assemble` method to behave as if the element were quadratic, creating "virtual" nodes on the midpoints of its edges and using 2nd-order shape functions. This is handled internally and does not require a different mesh file.
* **Public Functions**:
  * `~PhysicsField()`: Virtual destructor.
  * `setElementOrder(int order)` and `getElementOrder() const`: Sets and gets the mathematical order of the approximation (e.g., 1 for linear, 2 for quadratic).
  * `setup(Core::Mesh& mesh, Core::DOFManager& dof_manager)`: Sets up the field.
  * `assemble()`: Assembles the system matrices (K, M). This method is now much more powerful. It internally constructs a "virtual" higher-order element, gathers DOFs for both vertex and edge nodes, and uses the appropriate higher-order shape functions from `Utils::ShapeFunctions`.
  * `addBC(std::unique_ptr<Core::BoundaryCondition> bc)`: Adds a boundary condition.
  * `applyBCs()`: Applies all added boundary conditions. This method now intelligently consolidates BCs to ensure each degree of freedom is constrained only once.
  * ... (other functions remain the same) ...

### **Current1D/2D/3D**, **Heat1D/2D/3D**, **Magnetic1D/2D/3D**
* **Description**: These classes now all implement the advanced `assemble` method. They correctly gather DOFs for both vertex and higher-order edge nodes and use the `Utils::Quadrature` and `Utils::ShapeFunctions` libraries to perform accurate numerical integration for any specified element order.