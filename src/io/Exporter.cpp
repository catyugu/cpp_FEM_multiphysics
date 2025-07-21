#include "io/Exporter.hpp"
#include "core/Problem.hpp"
#include <core/mesh/Node.hpp>
#include <core/mesh/Element.hpp>
#include "utils/SimpleLogger.hpp"
#include <fstream>
#include <iomanip>
#include <core/DOFManager.hpp>
#include <physics/PhysicsField.hpp>

namespace IO {

bool Exporter::write_vtk(const std::string& filename, const Core::Problem& problem) {
    auto& logger = Utils::Logger::instance();
    logger.info("Exporting results to VTK file: ", filename);

    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        logger.error("Failed to open file for writing: ", filename);
        return false;
    }

    const auto& mesh = problem.getMesh();
    const auto& dof_manager = problem.getDofManager();
    const auto& fields = problem.getFields();

    // --- VTK Header ---
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "FEM Solver Results\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    // --- 1. POINTS Section: Write the geometric nodes ---
    const auto& nodes = mesh.getNodes();
    vtk_file << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        const auto& coords = node->getCoords();
        vtk_file << coords[0] << " " << coords[1] << " " << (coords.size() > 2 ? coords[2] : 0.0) << "\n";
    }
    vtk_file << "\n";

    // --- 2. CELLS Section: Define element connectivity using Node IDs ---
    const auto& elements = mesh.getElements();
    size_t cell_list_size = 0;
    for (const auto& elem : elements) {
        cell_list_size += (elem->getNodes().size() + 1);
    }
    vtk_file << "CELLS " << elements.size() << " " << cell_list_size << "\n";
    for (const auto& elem : elements) {
        const auto& elem_nodes = elem->getNodes();
        vtk_file << elem_nodes.size();
        for (const auto& node : elem_nodes) {
            vtk_file << " " << node->getId(); // Use the geometric node's ID
        }
        vtk_file << "\n";
    }
    vtk_file << "\n";

    // --- 3. CELL_TYPES Section: Define element geometry ---
    vtk_file << "CELL_TYPES " << elements.size() << "\n";
    for (const auto& elem : elements) {
        std::string typeName = elem->getTypeName();
        if (typeName == "LineElement") vtk_file << 3 << "\n";      // VTK_LINE
        else if (typeName == "TriElement") vtk_file << 5 << "\n";   // VTK_TRIANGLE
        else if (typeName == "TetElement") vtk_file << 10 << "\n";  // VTK_TETRA
        else vtk_file << 1 << "\n";                               // VTK_VERTEX
    }
    vtk_file << "\n";

    // --- 4. POINT_DATA Section: Write nodal solution values ---
    vtk_file << "POINT_DATA " << nodes.size() << "\n";
    if (fields.empty()) return true; // No data to write

    // Use the solution from any field, as the solver now makes them consistent.
    const auto& final_solution = fields[0]->getSolution();

    for (const auto& field : fields) {
        vtk_file << "SCALARS " << field->getVariableName() << " double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        std::string var_name = field->getVariableName();

        for (const auto& node : nodes) {
            int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
            if (dof_idx != -1 && dof_idx < final_solution.size()) {
                vtk_file << std::fixed << std::setprecision(6) << final_solution(dof_idx) << "\n";
            } else {
                vtk_file << 0.0 << "\n"; // Default value for nodes without this variable
            }
        }
    }

    vtk_file.close();
    logger.info("Successfully wrote VTK file: ", filename);
    return true;
}

} // namespace IO