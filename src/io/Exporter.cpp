#include "io/Exporter.hpp"
#include "core/Problem.hpp"
#include "core/Node.hpp"
#include "core/Element.hpp"
#include "utils/SimpleLogger.hpp"
#include <fstream>
#include <iomanip>
#include <core/DOFManager.hpp>
#include <core/Mesh.hpp>
#include <physics/PhysicsField.hpp>

namespace IO {

bool Exporter::write_vtk(const std::string& filename, const Core::Problem& problem)  {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Exporting results to VTK file: ", filename);

    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        logger.error("Failed to open file for writing: ", filename);
        return false;
    }

    const auto& mesh = problem.getMesh();
    const auto& dof_manager = problem.getDofManager();

    // 1. VTK Header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "FEM Solver Results\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    // 2. Points (Nodes)
    const auto& nodes = mesh.getNodes();
    vtk_file << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        const auto& coords = node->getCoords();
        vtk_file << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
    vtk_file << "\n";

    // 3. Cells (Elements)
    const auto& elements = mesh.getElements();
    size_t cell_list_size = 0;
    for (const auto& elem : elements) {
        cell_list_size += (elem->getNumNodes() + 1);
    }
    vtk_file << "CELLS " << elements.size() << " " << cell_list_size << "\n";
    for (const auto& elem : elements) {
        vtk_file << elem->getNumNodes();
        for (const auto& node : elem->getNodes()) {
            vtk_file << " " << node->getId();
        }
        vtk_file << "\n";
    }
    vtk_file << "\n";

    // 4. Cell Types
    vtk_file << "CELL_TYPES " << elements.size() << "\n";
    for (const auto& elem : elements) {
        // VTK_LINE has type ID 3
        if (std::string(elem->getTypeName()) == "LineElement") {
            vtk_file << 3 << "\n";
        } else {
            logger.warn("Unsupported element type for VTK export: ", elem->getTypeName());
            vtk_file << 1 << "\n"; // VTK_VERTEX as fallback
        }
    }
    vtk_file << "\n";

    // 5. Point Data (Nodal Results)
    vtk_file << "POINT_DATA " << nodes.size() << "\n";

    // Iterate through all solved physics fields and write their results
    for (const auto& var_name : dof_manager.getVariableNames()) {
        auto* field = problem.getField(var_name);
        if (field && field->getSolution().size() > 0) {
            vtk_file << "SCALARS " << var_name << " double 1\n";
            vtk_file << "LOOKUP_TABLE default\n";
            for (const auto& node : nodes) {
                int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
                if(dof_idx != -1) {
                    vtk_file << std::fixed << std::setprecision(6) << field->getSolution()(dof_idx) << "\n";
                } else {
                     vtk_file << 0.0 << "\n"; // Default value if DOF not found
                }
            }
        }
    }

    vtk_file.close();
    logger.info("Successfully wrote VTK file: ", filename);
    return true;
}

} // namespace IO
