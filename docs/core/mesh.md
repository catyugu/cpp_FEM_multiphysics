# **Core Mesh Components**

This document describes the fundamental components related to the mesh, including nodes and the various element types that form the geometric domain of the simulation.

## **Classes**

### **Mesh**

* **Description**: A container for all nodes and elements in the simulation. It also provides static factory methods for generating simple, uniform meshes.
* **Public Functions**:
    * \~Mesh()
    * addNode(Node\* node)
    * addElement(Element\* element)
    * getNode(int id) const: returns Node\*
    * getElement(int id) const: returns Element\*
    * getElements() const: returns const std::vector\<Element\*\>&
    * getNodes() const: returns const std::vector\<Node\*\>&
    * create\_uniform\_1d\_mesh(double length, int num\_elements): returns Mesh\*
    * create\_uniform\_2d\_mesh(double width, double height, int nx, int ny): returns Mesh\*
    * create\_uniform\_3d\_mesh(double width, double height, double depth, int nx, int ny, int nz): returns Mesh\*
* **Private Members**:
    * nodes\_: std::vector\<Node\*\>
    * elements\_: std::vector\<Element\*\>
    * node\_map\_: std::map\<int, Node\*\>
    * element\_map\_: std::map\<int, Element\*\>

### **Node**

* **Description**: Represents a point in space with an ID and coordinates.
* **Public Functions**:
    * Node(int id, double x, double y \= 0.0, double z \= 0.0)
    * getId() const: returns int
    * getCoords() const: returns const std::vector\<double\>&
* **Private Members**:
    * id\_: int
    * coords\_: std::vector\<double\>

### **Element**

* **Description**: An abstract base class for all mesh elements.
* **Public Functions**:
    * Element(int id)
    * getId() const: returns int
    * getNodes() const: returns const std::vector\<Node\*\>&
    * addNode(Node\* node)
    * getNumNodes() const: returns size\_t (pure virtual)
    * getTypeName() const: returns const char\* (pure virtual)
* **Protected Members**:
    * id\_: int
    * nodes\_: std::vector\<Node\*\>
* **Derived Classes**:
    * LineElement
    * TriElement
    * TetElement

### **LineElement**

* **Description**: A concrete implementation for a 1D line element with 2 nodes.
* **Public Functions**:
    * LineElement(int id)
    * getNumNodes() const override: returns size\_t
    * getTypeName() const override: returns const char\*
    * getLength() const: returns double

### **TriElement**

* **Description**: Represents a 2D, 3-node, linear triangular element.
* **Public Functions**:
    * TriElement(int id)
    * getNumNodes() const override: returns size\_t
    * getTypeName() const override: returns const char\*
    * getArea() const: returns double
    * getBMatrix() const: returns Eigen::Matrix\<double, 2, 3\>

### **TetElement**

* **Description**: Represents a 3D, 4-node, linear tetrahedral element.
* **Public Functions**:
    * TetElement(int id)
    * getNumNodes() const override: returns size\_t
    * getTypeName() const override: returns const char\*
    * getVolume() const: returns double
    * getBMatrix() const: returns Eigen::Matrix\<double, 3, 4\>