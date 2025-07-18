#include "io/Exporter.hpp"
#include "core/Problem.hpp"
#include <core/mesh/Node.hpp>
#include <core/mesh/Element.hpp>
#include "utils/SimpleLogger.hpp"
#include <fstream>
#include <iomanip>
#include <core/DOFManager.hpp>
#include <physics/PhysicsField.hpp>
#include <map>

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

    // Check if we have any higher-order fields. If so, we need a more complex export.
    bool has_higher_order_fields = false;
    if (!fields.empty()) {
        // We assume all fields use the same mesh and DOF manager, so we can use the first field to check.
        const auto& first_field = fields[0];
        if (first_field->getElementOrder() > 1) {
            has_higher_order_fields = true;
        }
    }

    // --- Header ---
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "FEM Solver Results\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    // --- Points Section ---
    // For higher-order, the "points" in the VTK file correspond to the DOFs, not just mesh nodes.
    // We must generate coordinates for all DOFs.
    std::map<int, std::vector<double>> dof_coords;
    for (const auto& var_name : dof_manager.getVariableNames()) {
        // 1. Get coordinates for vertex DOFs
        for (const auto& node : mesh.getNodes()) {
            int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
            if (dof_idx != -1) {
                dof_coords[dof_idx] = node->getCoords();
            }
        }
        // 2. Generate coordinates for edge DOFs
        if (has_higher_order_fields) {
            for (const auto& elem : mesh.getElements()) {
                 const auto& elem_nodes = elem->getNodes();
                 std::vector<std::pair<int, int>> edges;
                 if (elem->getDimension() == 1 && elem_nodes.size() == 2) edges = {{0, 1}};
                 if (elem->getDimension() == 2 && elem_nodes.size() == 3) edges = {{0, 1}, {1, 2}, {2, 0}};
                 if (elem->getDimension() == 3 && elem_nodes.size() == 4) edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};

                 for (const auto& edge : edges) {
                    auto* n1 = elem_nodes[edge.first];
                    auto* n2 = elem_nodes[edge.second];
                    std::vector<int> edge_node_ids = {n1->getId(), n2->getId()};
                    int dof_idx = dof_manager.getEdgeEquationIndex(edge_node_ids, var_name);
                    if (dof_idx != -1 && dof_coords.find(dof_idx) == dof_coords.end()) {
                        std::vector<double> mid_coords = {
                            (n1->getCoords()[0] + n2->getCoords()[0]) / 2.0,
                            (n1->getCoords()[1] + n2->getCoords()[1]) / 2.0,
                            (n1->getCoords()[2] + n2->getCoords()[2]) / 2.0
                        };
                        dof_coords[dof_idx] = mid_coords;
                    }
                }
            }
        }
    }

    vtk_file << "POINTS " << dof_manager.getNumEquations() << " double\n";
    for (size_t i = 0; i < dof_manager.getNumEquations(); ++i) {
        if (dof_coords.count(i)) {
             vtk_file << dof_coords[i][0] << " " << dof_coords[i][1] << " " << dof_coords[i][2] << "\n";
        } else {
             // This case should ideally not happen if all DOFs are accounted for.
             vtk_file << "0.0 0.0 0.0\n";
        }
    }
    vtk_file << "\n";


    // --- Cells Section ---
    const auto& elements = mesh.getElements();
    size_t cell_list_size = 0;
    for (const auto& elem : elements) {
        // Use the first field to get element info, assuming all fields on an element have the same order.
        if (!fields.empty()) {
            elem->setOrder(fields[0]->getElementOrder());
        }
        cell_list_size += (elem->getNumNodes() + 1);
    }
    vtk_file << "CELLS " << elements.size() << " " << cell_list_size << "\n";
    for (const auto& elem : elements) {
         if (fields.empty()) continue;
         const auto& field = fields[0]; // Use first field to get DOFs
         elem->setOrder(field->getElementOrder());
         auto dofs = field->get_element_dofs(elem);
         vtk_file << dofs.size();
         for(int dof_idx : dofs) {
             vtk_file << " " << dof_idx;
         }
         vtk_file << "\n";
    }
    vtk_file << "\n";

    // --- Cell Types Section (FIXED) ---
    vtk_file << "CELL_TYPES " << elements.size() << "\n";
    for (const auto& elem : elements) {
        int order = fields.empty() ? 1 : fields[0]->getElementOrder();
        std::string typeName = elem->getTypeName();
        if (typeName == "LineElement") {
            vtk_file << (order > 1 ? 21 : 3) << "\n"; // VTK_QUADRATIC_LINE or VTK_LINE
        } else if (typeName == "TriElement") {
            vtk_file << (order > 1 ? 22 : 5) << "\n"; // VTK_QUADRATIC_TRIANGLE or VTK_TRIANGLE
        } else if (typeName == "TetElement") {
            vtk_file << (order > 1 ? 24 : 10) << "\n"; // VTK_QUADRATIC_TETRA or VTK_TETRA
        } else {
            logger.warn("Unsupported element type for VTK export: ", typeName);
            vtk_file << 1 << "\n"; // VTK_VERTEX as fallback
        }
    }
    vtk_file << "\n";

    // --- Point Data Section (FIXED) ---
    vtk_file << "POINT_DATA " << dof_manager.getNumEquations() << "\n";
    for (const auto& field : fields) {
        if (field->getSolution().size() > 0) {
            vtk_file << "SCALARS " << field->getVariableName() << " double 1\n";
            vtk_file << "LOOKUP_TABLE default\n";
            // Simply write out the entire solution vector, as points now map directly to DOFs
            for (int i = 0; i < field->getSolution().size(); ++i) {
                vtk_file << std::fixed << std::setprecision(6) << field->getSolution()(i) << "\n";
            }
        }
    }

    vtk_file.close();
    logger.info("Successfully wrote VTK file: ", filename);
    return true;
}

} // namespace IO