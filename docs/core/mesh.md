# **Core Mesh Components**

This document describes the fundamental components related to the mesh, including nodes, the various element types that form the geometric domain of the simulation, and the new **Finite Element Values (FEValues)** calculator which is central to the p-refinement strategy.

---
## **Classes**

### **Mesh**

* **Description**: A container for all nodes and elements in the simulation. It also provides static factory methods for generating simple, uniform meshes.
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

### **ElementGeometry**

* **Description**: A new class that holds the pure geometric information of an element's vertices. It is a fundamental component for the `FEValues` calculator.
* **Public Functions**:
  * `ElementGeometry(const std::vector<Node*>& vertex_nodes, int dimension)`: Constructor that takes a list of vertex nodes and the element's dimension.
  * `get_vertex_coords() const`: Returns an Eigen matrix of the vertex coordinates.
  * `get_dimension() const`: Returns the intrinsic dimension of the element (1, 2, or 3).
  * `get_num_vertices() const`: Returns the number of geometric vertices.
* **Private Members**:
  * `vertex_nodes_`: `std::vector<Node*>`
  * `vertex_coords_`: `Eigen::MatrixXd`
  * `dimension_`: `int`

### **FEValues**

* **Description**: The **Finite Element Values calculator**. This crucial new class is responsible for pre-calculating and caching all values needed for FEM assembly at each quadrature point. This includes shape function values (`N`), their gradients in real coordinates (`∇N`), and the Jacobian determinant multiplied by the quadrature weight (`detJ * w_q`).
* **Key Concept**: By pre-calculating these values, `FEValues` significantly simplifies the `assemble` methods within the `PhysicsField` classes, making the assembly loop incredibly fast and clean.
* **Public Functions**:
  * `FEValues(const ElementGeometry& geom, int order, int quad_order)`: Constructor that takes an `ElementGeometry` object, the desired mathematical `order` of the shape functions, and the `quad_order` for numerical integration.
  * `reinit(int q_point_index)`: Re-initializes the object to provide values for a specific quadrature point, making the cached values available via accessors.
  * `num_quadrature_points() const`: Returns the total number of quadrature points for this element and quadrature order.
  * `get_shape_values() const`: Returns the shape function values (`N`) at the current quadrature point.
  * `get_shape_gradients() const`: Returns the shape function gradients in the real coordinate system (`∇N`), often referred to as the B-matrix, at the current quadrature point.
  * `get_detJ_times_weight() const`: Returns the product of the Jacobian determinant and the quadrature weight at the current quadrature point.
* **Private Members**:
  * `N_values_`: `Eigen::VectorXd`
  * `dN_dx_values_`: `Eigen::MatrixXd`
  * `detJ_x_w_`: `double`
  * `geometry_`: `const ElementGeometry&`
  * `fe_order_`: `int`
  * `quadrature_points_`: `std::vector<Utils::QuadraturePoint>`
  * `all_N_values_`: `std::vector<Eigen::VectorXd>`
  * `all_dN_dx_values_`: `std::vector<Eigen::MatrixXd>`
  * `all_detJ_x_w_`: `std::vector<double>`

### **Element**

* **Description**: An abstract base class for all mesh elements. The mesh itself is always composed of geometrically linear elements (e.g., a triangle always has 3 vertex nodes). This class has been redesigned to integrate with `ElementGeometry` and `FEValues`, becoming the primary source of finite element calculation data.
* **Public Functions**:
  * `Element(int id)`
  * `getId() const`: returns `int`
  * `getNodes() const`: returns `const std::vector<Node*>&`
  * `addNode(Node* node)`
  * `getNumNodes() const`: Returns the number of nodes based on the element's **mathematical order** (`order_`), not just its geometric vertices.
  * `getTypeName() const`: returns `const char*` (pure virtual)
  * `getDimension() const`: Pure virtual function that returns the intrinsic physical dimension of the element (1, 2, or 3).
  * `setOrder(int order)` and `getOrder() const`: Set and get the mathematical order for the element's approximation.
  * `update_geometry()`: Helper method to initialize or update the internal `ElementGeometry` object once nodes have been added.
  * `create_fe_values(int quad_order)`: The primary new role of `Element`. It creates and returns a `std::unique_ptr<FEValues>` object for this element, using its own `order_` and the specified `quad_order`.
* **Protected Members**:
  * `id_`: `int`
  * `nodes_`: `std::vector<Node*>`
  * `order_`: `int` (Defaults to 1 for linear)
  * `geometry_`: `std::unique_ptr<ElementGeometry>` (Each element now holds its own geometry).
* **Derived Classes**:
  * `LineElement`
  * `TriElement`
  * `TetElement`

### **LineElement**, **TriElement**, **TetElement**

* **Description**: Concrete implementations for linear mesh elements. Their `getNumNodes()` method now dynamically returns the correct number of nodes required for higher-order formulations (e.g., 10 for a quadratic tetrahedron) based on the `order_` property set on the base `Element` class.