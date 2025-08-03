# **Core Mesh Components**

This document describes the fundamental components related to the mesh, including nodes, the various element types that form the geometric domain of the simulation, and the **Reference Element Cache** which is central to the p-refinement strategy.

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

* **Description**: A class that holds the pure geometric information of an element's vertices. It is a fundamental component for the `FEValues` calculator.
* **Public Functions**:
  * `ElementGeometry(const std::vector<Node*>& vertex_nodes, int dimension)`: Constructor that takes a list of vertex nodes and the element's dimension.
  * `get_vertex_coords() const`: Returns an Eigen matrix of the vertex coordinates.
  * `get_dimension() const`: Returns the intrinsic dimension of the element (1, 2, or 3).
  * `get_num_vertices() const`: Returns the number of geometric vertices.
* **Private Members**:
  * `vertex_nodes_`: `std::vector<Node*>`
  * `vertex_coords_`: `Eigen::MatrixXd`
  * `dimension_`: `int`

### **ReferenceElementCache** and **FEValues**

* **Description**: This pair of classes represents a significant architectural improvement. Instead of each `Element` creating its own `FEValues` object from scratch, a static, thread-safe `ReferenceElementCache` is used to pre-compute and store geometry-independent data.
* **`ReferenceElementCache`**:
  * A singleton-like class that stores `ReferenceElementData` (quadrature points, shape function values, and their derivatives in natural coordinates).
  * It uses a cache key of `{element_type, fe_order, quad_order}` to ensure that this data is calculated only once per simulation run.
* **`FEValues`**:
  * A lightweight "calculator" or "view" that is instantiated inside the `assemble` loop.
  * Its constructor takes the specific `ElementGeometry` and a reference to the pre-computed `ReferenceElementData` from the cache.
  * It then performs the geometry-dependent calculations (Jacobians, gradients in real coordinates) for that specific element.
  * This design dramatically speeds up the assembly process by avoiding redundant calculations of shape functions and quadrature rules for elements of the same type and order.
* **Public Functions (`FEValues`)**:
  * `FEValues(const ElementGeometry& geom, int order, const ReferenceElementData& ref_data)`: The lightweight constructor.
  * `reinit(int q_point_index)`: Sets the calculator to a specific quadrature point.
  * `num_quadrature_points() const`: Returns the number of quadrature points.
  * `get_shape_values() const`: Returns the shape function values (`N`).
  * `get_shape_gradients() const`: Returns the shape function gradients in the real coordinate system (`∇N`).
  * `get_detJ_times_weight() const`: Returns `det(J) * w_q`.

### **Element**

* **Description**: An abstract base class for all mesh elements. The mesh is always composed of geometrically linear elements.
* **Public Functions**:
  * `Element(int id)`
  * `getId() const`
  * `getNodes() const`: Returns the geometric vertex nodes of the element.
  * `getNumNodes() const`: Returns the number of nodes based on the element's **mathematical order** (`order_`), not just its geometric vertices (e.g., a quadratic tet has 10 nodes).
  * `getTypeName() const`: Returns a string like "TriElement" or "TetElement", used as a key for the `ReferenceElementCache`.
  * `getDimension() const`: Pure virtual function returning the element's dimension (1, 2, or 3).
  * `setOrder(int order)` and `getOrder() const`: Set and get the mathematical order for the element's approximation.
  * `getGeometry()`: Returns a reference to the `ElementGeometry` object.
  * `createFEValues(int quad_order)`: Returns a `unique_ptr<FEValues>` ptr for the element.
* **Protected Members**:
  * `id_`: `int`
  * `nodes_`: The geometric vertices (`std::vector<Node*>`).
  * `order_`: The mathematical approximation order (`int`).
  * `geometry_`: `std::unique_ptr<ElementGeometry>`.

### **LineElement**, **TriElement**, **TetElement**

* **Description**: Concrete implementations for linear mesh elements. Their `getNumNodes()` method dynamically returns the correct number of nodes required for higher-order formulations based on the `order_` property. The `TetElement::getVolume()` method now correctly returns a **signed volume**, which is used by the importer to detect and fix inverted elements.