# **Utils Namespace**

This namespace provides general-purpose utilities that support the core numerical methods of the simulation framework.

---
## **Classes**

### **ShapeFunctions**
* **Description**: A static utility class that provides the basis functions (and their derivatives) required for finite element formulations of different orders (p-refinement). It provides the mathematical foundation for approximating the solution over an element.
* **Key Static Functions**:
    * `getLineShapeFunctions(int order, ...)`
    * `getTriShapeFunctions(int order, ...)`
    * `getTetShapeFunctions(int order, ...)`
    * `getLineShapeFunctionDerivatives(int order, ...)`
    * `getTriShapeFunctionDerivatives(int order, ...)`
    * `getTetShapeFunctionDerivatives(int order, ...)`
* **Supported Orders**: The class currently provides shape functions for linear and quadratic (order 1 and 2) elements for lines, triangles, and tetrahedra.

### **Quadrature**
* **Description**: A static utility class that provides sets of points and corresponding weights for performing accurate numerical integration (Gaussian Quadrature) over elements. The choice of quadrature order is critical for correctly integrating the element matrices, especially for higher-order elements.
* **Key Static Functions**:
    * `getLineQuadrature(int order)`
    * `getTriangleQuadrature(int order)`
    * `getTetrahedronQuadrature(int order)`
* **Functionality**: Returns a `std::vector` of `QuadraturePoint` structs, each containing the coordinate of the integration point in the element's natural coordinate system and its corresponding weight.

### **SimpleLogger**
* **Description**: A thread-safe, singleton logger for writing timestamped and color-coded messages to the console and/or a log file.

### **Exceptions**
* **Description**: A namespace containing custom exception types (`FileIOException`, `SolverException`, `ConfigurationException`) to provide more specific error information than standard exceptions.