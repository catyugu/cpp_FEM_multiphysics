# **Utils Namespace**

This namespace provides general-purpose utilities that support the core numerical methods of the simulation framework.

---
## **Classes**

### **ShapeFunctions**
* **Description**: A static utility class that provides the basis functions (and their derivatives) required for finite element formulations of different orders (p-refinement). It provides the mathematical foundation for approximating the solution over an element.
* **Key Static Functions**:
  * `getLineShapeFunctions(int order, double xi)`: Returns shape function values for 1D line elements.
  * `getTriShapeFunctions(int order, double xi, double eta)`: Returns shape function values for 2D triangle elements.
  * `getTetShapeFunctions(int order, double xi, double eta, double zeta)`: Returns shape function values for 3D tetrahedron elements.
  * `getLineShapeFunctionDerivatives(int order, double xi)`: Returns shape function derivatives for 1D line elements.
  * `getTriShapeFunctionDerivatives(int order, double xi, double eta)`: Returns shape function derivatives for 2D triangle elements. Includes a **fix** for the derivative of the N5 shape function (`4*L2*L3`) with respect to `eta` for order 2.
  * `getTetShapeFunctionDerivatives(int order, double xi, double eta, double zeta)`: Returns shape function derivatives for 3D tetrahedron elements.
* **Supported Orders**: The class currently provides shape functions for linear (order 1) and quadratic (order 2) elements for lines, triangles, and tetrahedra. It supports up to order 5 for line elements.

### **Quadrature**
* **Description**: A static utility class that provides sets of points and corresponding weights for performing accurate numerical integration (Gaussian Quadrature) over elements. The choice of quadrature order is critical for correctly integrating the element matrices, especially for higher-order elements.
* **Key Static Functions**:
  * `getLineQuadrature(int order)`: Returns quadrature points and weights for 1D line elements. Now supports up to order 5.
  * `getTriangleQuadrature(int order)`: Returns quadrature points and weights for 2D triangle elements.
  * `getTetrahedronQuadrature(int order)`: Returns quadrature points and weights for 3D tetrahedron elements.
* **Functionality**: Returns a `std::vector` of `QuadraturePoint` structs, each containing the coordinate of the integration point in the element's natural coordinate system and its corresponding weight.

### **SimpleLogger**
* **Description**: A thread-safe, singleton logger for writing timestamped and color-coded messages to the console and/or a log file.

### **Exceptions**
* **Description**: A namespace containing custom exception types (`FileIOException`, `SolverException`, `ConfigurationException`) to provide more specific error information than standard exceptions.