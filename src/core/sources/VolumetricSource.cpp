#include "core/sources/VolumetricSource.hpp"
#include "core/DOFManager.hpp"
#include "core/mesh/Mesh.hpp"
#include "core/mesh/Element.hpp"
#include "core/mesh/LineElement.hpp"
#include "core/mesh/TriElement.hpp"
#include "core/mesh/TetElement.hpp"
#include <vector>
#include "core/FEValues.hpp" // Include FEValues

namespace Core {

    VolumetricSource::VolumetricSource(int element_id, double total_power, const std::string& tag)
        : SourceTerm(tag), element_id_(element_id), total_power_(total_power) {}

    void VolumetricSource::apply(Eigen::VectorXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const {
        Element* elem = mesh.getElement(element_id_);
        if (!elem) return;

        elem->setOrder(element_order);

        // 使用Element的缓存机制获取FEValues对象
        FEValues* fe_values = elem->getFEValues(element_order, AnalysisType::SCALAR_DIFFUSION);
        
        std::vector<int> element_dofs;
        const auto& vertex_nodes_of_element = elem->getNodes();
        const size_t num_vertices = vertex_nodes_of_element.size();

        element_dofs.reserve(elem->getNumNodes()); // Pre-allocate to prevent reallocations

        // 1. Get Vertex DOFs
        for (size_t i = 0; i < num_vertices; ++i) {
            element_dofs.push_back(dof_manager.getEquationIndex(vertex_nodes_of_element[i]->getId(), var_name));
        }

        // 2. Get Higher-Order DOFs (if any) in their correct canonical order
        if (element_order > 1) {
            if (auto* line_elem = dynamic_cast<const Core::LineElement*>(elem)) {
                element_dofs.insert(element_dofs.begin() + 1, dof_manager.getEdgeEquationIndex({line_elem->getNodes()[0]->getId(), line_elem->getNodes()[1]->getId()}, var_name));
            }
            else if (auto* tri_elem = dynamic_cast<const Core::TriElement*>(elem)) {
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto& edge : edges) {
                    element_dofs.push_back(dof_manager.getEdgeEquationIndex({tri_elem->getNodes()[edge.first]->getId(), tri_elem->getNodes()[edge.second]->getId()}, var_name));
                }
            } else if (auto* tet_elem = dynamic_cast<const Core::TetElement*>(elem)) {
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                for (const auto& edge : edges) {
                    element_dofs.push_back(dof_manager.getEdgeEquationIndex({tet_elem->getNodes()[edge.first]->getId(), tet_elem->getNodes()[edge.second]->getId()}, var_name));
                }
            }
        }
        // This 'element_dofs' now contains all DOFs for the element, correctly ordered according to ShapeFunctions convention.


        // Calculate the integral of each shape function over the element volume/area
        // Sum of integral_Ni_dV should equal the element's volume/area
        std::vector<double> integral_Ni_dV(elem->getNumNodes(), 0.0);
        double total_element_volume_recomputed = 0.0;

        for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
            fe_values->reinit(q_p);
            const auto& N = fe_values->get_shape_values();
            const double detJ_x_w = fe_values->get_detJ_times_weight();

            for (size_t i = 0; i < N.size(); ++i) {
                integral_Ni_dV[i] += N(i) * detJ_x_w;
            }
            total_element_volume_recomputed += detJ_x_w; // Sum of detJ*w is the element volume/area
        }

        // Distribute total_power_ proportionally based on the integral of shape functions
        if (std::abs(total_element_volume_recomputed) < 1e-12) { // Avoid division by zero for degenerate elements
            // Utils::Logger::instance().warn("Element ", element_id_, " has zero or near-zero volume, skipping source application.");
            return;
        }

        for (size_t i = 0; i < element_dofs.size(); ++i) {
            if (element_dofs[i] != -1) { // Ensure the DOF actually exists and is mapped
                // The nodal force contribution is P_total * (integral_Ni_dV / Element_Volume)
                F(element_dofs[i]) += total_power_ * (integral_Ni_dV[i] / total_element_volume_recomputed);
            }
        }
    }
} // namespace Core