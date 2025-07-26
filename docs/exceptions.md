# **Exceptions Namespace**

This document describes the custom exception classes used in the C++ FEM Multiphysics project for robust error handling.

---
## **Exception Classes**

### **FileIOException**
* **Description**: Thrown when a file I/O error occurs, such as failing to open a mesh file or a results file.
* **Inherits from**: `std::runtime_error`

### **SolverException**
* **Description**: Thrown when a numerical solver fails. This can happen if a direct solver encounters a singular matrix, or if an iterative solver fails to converge within the specified number of iterations.
* **Inherits from**: `std::runtime_error`

### **ConfigurationException**
* **Description**: Thrown when there is an error in the problem setup. Examples include requesting a material property that hasn't been defined, using an unsupported element type, or trying to solve a problem before calling the `setup()` method.
* **Inherits from**: `std::runtime_error`