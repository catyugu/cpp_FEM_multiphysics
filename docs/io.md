# **IO Namespace**

This namespace contains the specific implementations for different physical simulations.

---
## **Classes**

### **Exporter**
* **Description**: A utility class for exporting simulation results to various file formats.
* **Public Functions**:
    * `write_vtk(const std::string& filename, const Core::Problem& problem)`: Writes the mesh and all solved nodal data from a Problem to a legacy VTK file.

### **Importer**
* **Description**: A utility class for importing mesh data from various file formats.
* **Public Functions**:
    * `read_comsol_mphtxt(const std::string& filename)`: Reads a mesh from a COMSOL .mphtxt text file.