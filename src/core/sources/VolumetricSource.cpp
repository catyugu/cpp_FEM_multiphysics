#include "core/sources/VolumetricSource.hpp"
#include "core/DOFManager.hpp"
#include "core/mesh/Mesh.hpp"
#include "core/mesh/Element.hpp"
#include "core/mesh/LineElement.hpp" // Required for dynamic_cast
#include "core/mesh/TriElement.hpp"   // Required for dynamic_cast
#include "core/mesh/TetElement.hpp"   // Required for dynamic_cast
#include <vector> // Required for std::vector

namespace Core {

    VolumetricSource::VolumetricSource(int element_id, double total_power, const std::string& tag)
        : SourceTerm(tag), element_id_(element_id), total_power_(total_power) {}

    // MODIFICATION: Added 'int element_order' to the signature
    void VolumetricSource::apply(Eigen::MatrixXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const {
        Element* elem = mesh.getElement(element_id_);
        if (!elem) return;

        elem->setOrder(element_order); // Ensure element knows its current order for getNumNodes()

        std::vector<int> element_dofs;
        const auto& vertex_nodes_of_element = elem->getNodes();
        const size_t num_vertices = vertex_nodes_of_element.size();

        // 1. Get Vertex DOFs
        for (size_t i = 0; i < num_vertices; ++i) {
            element_dofs.push_back(dof_manager.getEquationIndex(vertex_nodes_of_element[i]->getId(), var_name));
        }

        // 2. Get Higher-Order DOFs (if any)
        if (element_order > 1) {
            if (auto* line_elem = dynamic_cast<const Core::LineElement*>(elem)) {
                // Canonical order for P2 Line: v0, midpoint, v1
                element_dofs.push_back(dof_manager.getEdgeEquationIndex({line_elem->getNodes()[0]->getId(), line_elem->getNodes()[1]->getId()}, var_name));
            }
            else if (auto* tri_elem = dynamic_cast<const Core::TriElement*>(elem)) {
                // Canonical edge order for Tri6: edge(0,1), edge(1,2), edge(2,0)
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto& edge : edges) {
                    element_dofs.push_back(dof_manager.getEdgeEquationIndex({tri_elem->getNodes()[edge.first]->getId(), tri_elem->getNodes()[edge.second]->getId()}, var_name));
                }
            } else if (auto* tet_elem = dynamic_cast<const Core::TetElement*>(elem)) {
                // Canonical edge order for Tet10: edge(0,1), edge(0,2), edge(0,3), edge(1,2), edge(1,3), edge(2,3)
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                for (const auto& edge : edges) {
                    element_dofs.push_back(dof_manager.getEdgeEquationIndex({tet_elem->getNodes()[edge.first]->getId(), tet_elem->getNodes()[edge.second]->getId()}, var_name));
                }
            }
        }
        // This 'element_dofs' now contains all DOFs for the element, correctly ordered according to ShapeFunctions convention.

        // Distribute total_power_ evenly among all relevant DOFs for this element.
        // The total number of DOFs for the element is elem->getNumNodes() for the given element_order.
        // It's crucial that element_dofs.size() matches elem->getNumNodes() for consistency.
        if (element_dofs.size() == 0) return; // Avoid division by zero if no DOFs found (shouldn't happen)

        double nodal_share_per_dof = total_power_ / element_dofs.size();

        for (int dof_index : element_dofs) {
            if (dof_index != -1) { // Ensure the DOF actually exists and is mapped
                F(dof_index) += nodal_share_per_dof;
            }
        }
    }

} // namespace Core