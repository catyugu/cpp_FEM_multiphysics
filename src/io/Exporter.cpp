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

    // Header and Points sections remain the same...
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "FEM Solver Results\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";
    const auto& nodes = mesh.getNodes();
    vtk_file << "POINTS " << nodes.size() << " double\n";
    for (const auto& node : nodes) {
        const auto& coords = node->getCoords();
        vtk_file << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
    vtk_file << "\n";

    // Cells section remains the same...
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

    // --- FIX: Add support for TetElement type ---
    vtk_file << "CELL_TYPES " << elements.size() << "\n";
    for (const auto& elem : elements) {
        std::string typeName = elem->getTypeName();
        if (typeName == "LineElement") {
            vtk_file << 3 << "\n"; // VTK_LINE
        } else if (typeName == "TriElement") {
            vtk_file << 5 << "\n"; // VTK_TRIANGLE
        } else if (typeName == "TetElement") {
            vtk_file << 10 << "\n"; // VTK_TETRA
        } else {
            logger.warn("Unsupported element type for VTK export: ", typeName);
            vtk_file << 1 << "\n"; // VTK_VERTEX as fallback
        }
    }
    vtk_file << "\n";

    // Point Data section remains the same...
    vtk_file << "POINT_DATA " << nodes.size() << "\n";
    for (const auto& var_name : dof_manager.getVariableNames()) {
        const auto* field = problem.getField(var_name);
        if (field && field->getSolution().size() > 0) {
            vtk_file << "SCALARS " << var_name << " double 1\n";
            vtk_file << "LOOKUP_TABLE default\n";
            for (const auto& node : nodes) {
                int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
                if(dof_idx != -1) {
                    vtk_file << std::fixed << std::setprecision(6) << field->getSolution()(dof_idx) << "\n";
                } else {
                     vtk_file << 0.0 << "\n";
                }
            }
        }
    }

    vtk_file.close();
    logger.info("Successfully wrote VTK file: ", filename);
    return true;
}

} // namespace IO
