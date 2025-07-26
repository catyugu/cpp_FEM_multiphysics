# **IO Namespace**

This namespace provides utilities for importing and exporting data, primarily handling mesh files and simulation results. It is designed to be flexible, supporting common formats from meshing software like Gmsh and COMSOL, and visualization tools like ParaView.

---
## **Classes**

### **Importer**
* **Description**: A utility class containing static methods for importing mesh data and validation results from various file formats.
* **Public Static Functions**:
  * `read_comsol_mphtxt(const std::string& filename)`: Reads a mesh from a COMSOL `.mphtxt` text file.
    * **Functionality**: Intelligently detects the spatial dimension (2D or 3D) from the file and only reads the corresponding element types (triangles for 2D, tetrahedra for 3D). For 3D meshes, it automatically checks the signed volume of each tetrahedron and reorders nodes if necessary to fix inverted elements, ensuring a valid mesh for simulation.
  * `read_gmsh_msh(const std::string& filename)`: Reads a mesh from a Gmsh `.msh` file (ASCII format v2).
    * **Functionality**: Parses a `.msh` file to read nodes and elements (lines, triangles, tetrahedra), creating a `Core::Mesh` object.
  * `read_vtu_data(const std::string& filename, const std::string& data_array_name)`: Reads a single data array from a VTK Unstructured Grid (`.vtu`) file.
    * **Functionality**: Primarily used in testing to load a single vector of reference data for validation.
  * `read_vtu_points_and_data(const std::string& filename, const std::vector<std::string>& data_array_names)`: Reads node coordinates and multiple specified data arrays from a `.vtu` file.
    * **Functionality**: Returns a `VtuData` struct containing point coordinates and a map of the requested data arrays. This is the primary tool for validation against reference results from other software like COMSOL.

---
### **Exporter**
* **Description**: A utility class for exporting simulation results to a format suitable for post-processing.
* **Public Static Functions**:
  * `write_vtk(const std::string& filename, const Core::Problem& problem)`: Writes the mesh and all solved nodal data from a `Problem` to a legacy VTK (`.vtk`) file.
    * **Functionality**: Exports the simulation's mesh and the solution data for each physics field. It correctly handles different element types (Line, Triangle, Tetrahedron) by mapping them to their corresponding VTK cell types. It writes both **scalar** and **vector** nodal (point) data, allowing for visualization of fields like Temperature, Voltage, and Magnetic Vector Potential in tools like ParaView or VisIt.