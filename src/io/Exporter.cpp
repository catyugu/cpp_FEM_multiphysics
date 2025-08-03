// src/io/Exporter.cpp (Upgraded Version)

#include "io/Exporter.hpp"
#include "core/Problem.hpp"
#include <core/mesh/Node.hpp>
#include <core/mesh/Element.hpp>
#include "utils/SimpleLogger.hpp"
#include <fstream>
#include <iomanip>
#include <core/DOFManager.hpp>
#include <physics/PhysicsField.hpp>
#include "post/PostProcessor.hpp" // Ensure this is included

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
    const auto& post_results = problem.getPostProcessingResults(); // Get post-processing results

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
            vtk_file << " " << node->getId();
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
    if (!fields.empty()) {
        for (const auto& field : fields) {
            if (!field->isEnabled()) continue;
            std::string var_name = field->getVariableName();
            int num_components = field->getNumComponents();
            const auto& solution = field->getSolution();

            if (num_components == 1) {
                vtk_file << "SCALARS " << var_name << " double 1\n";
                vtk_file << "LOOKUP_TABLE default\n";
                for (const auto& node : nodes) {
                    int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
                    if (dof_idx != -1 && dof_idx < solution.size()) {
                        vtk_file << std::fixed << std::setprecision(6) << solution(dof_idx) << "\n";
                    } else {
                        vtk_file << 0.0 << "\n";
                    }
                }
            } else if (num_components > 1) { // Handle vectors (e.g., MagneticVectorPotential)
                vtk_file << "VECTORS " << var_name << " double\n";
                for (const auto& node : nodes) {
                    int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
                    if (dof_idx != -1 && (dof_idx + num_components - 1) < solution.size()) {
                        for(int c = 0; c < num_components; ++c) {
                            vtk_file << std::fixed << std::setprecision(6) << solution(dof_idx + c) << (c == num_components - 1 ? "" : " ");
                        }
                        if (num_components < 3) { // VTK VECTORS must have 3 components
                           for (int c = num_components; c < 3; ++c) vtk_file << " 0.0";
                        }
                        vtk_file << "\n";
                    } else {
                        vtk_file << "0.0 0.0 0.0\n";
                    }
                }
            }
        }
    }

    // --- 5. CELL_DATA Section: Write post-processed element values ---
    if (!post_results.empty()) {
        vtk_file << "CELL_DATA " << elements.size() << "\n";
        for (const auto& pair : post_results) {
            const auto& result = pair.second;

            if (result.dimension == 1) { // Scalar cell data
                vtk_file << "SCALARS " << result.name << " double 1\n";
                vtk_file << "LOOKUP_TABLE default\n";
            } else { // Vector cell data
                vtk_file << "VECTORS " << result.name << " double\n";
            }

            for (size_t i = 0; i < result.data.size(); ++i) {
                Eigen::VectorXd avg_val = Eigen::VectorXd::Zero(result.dimension);
                if (!result.data[i].empty()) {
                    for (const auto& qp_val : result.data[i]) {
                        avg_val += qp_val;
                    }
                    avg_val /= result.data[i].size();
                }

                // Write the averaged value for the element
                for (int d = 0; d < result.dimension; ++d) {
                    vtk_file << std::fixed << std::setprecision(6) << avg_val(d) << (d == result.dimension - 1 ? "" : " ");
                }
                if (result.dimension < 3 && result.dimension > 1) { // Pad vectors to 3 components for VTK
                    for (int d = result.dimension; d < 3; ++d) vtk_file << " 0.0";
                }
                vtk_file << "\n";
            }
        }
    }


    vtk_file.close();
    logger.info("Successfully wrote VTK file: ", filename);
    return true;
}

} // namespace IO