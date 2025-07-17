# **Exceptions Namespace**

This document describes the custom exception classes used in the C++ FEM Multiphysics project.

---
## **Exception Classes**

### **FileIOException**
* **Description**: Thrown when a file I/O error occurs, such as failing to open a mesh file.
* **Inherits from**: `std::runtime_error`

### **SolverException**
* **Description**: Thrown when a numerical solver fails, for example, due to a singular matrix, or if an iterative solver does not converge.
* **Inherits from**: `std::runtime_error`

### **ConfigurationException**
* **Description**: Thrown when there is an error in the problem configuration, such as a missing material property, an unsupported element type for FEValues, or an uninitialized solver.
* **Inherits from**: `std::runtime_error`