# **Core Mesh Components**

This document describes the fundamental components related to the mesh, including nodes and the various element types that form the geometric domain of the simulation.

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

### **Element**

* **Description**: An abstract base class for all mesh elements. Note that the mesh itself is always composed of geometrically linear elements (e.g., a triangle always has 3 vertex nodes).
* **Public Functions**:
  * `Element(int id)`
  * `getId() const`: returns `int`
  * `getNodes() const`: returns `const std::vector<Node*>&`
  * `addNode(Node* node)`
  * `getNumNodes() const`: Returns the number of nodes based on the element's mathematical order.
  * `getTypeName() const`: returns `const char*` (pure virtual)
  * `setOrder(int order)` and `getOrder() const`: Get and set the mathematical order for the element's formulation.
* **Protected Members**:
  * `id_`: `int`
  * `nodes_`: `std::vector<Node*>`
  * `order_`: `int` (Defaults to 1 for linear)
* **Derived Classes**:
  * `LineElement`
  * `TriElement`
  * `TetElement`

### **LineElement**, **TriElement**, **TetElement**

* **Description**: Concrete implementations for linear mesh elements. Their `getNumNodes()` method is now capable of returning the correct number of nodes for higher-order formulations (e.g., 10 for a quadratic tetrahedron) based on the `order_` property.