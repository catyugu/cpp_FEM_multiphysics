#include "core/sources/PrescribedCurrentDensity.hpp"
#include "core/mesh/Element.hpp"
#include "core/FEValues.hpp"
#include "physics/PhysicsField.hpp"
#include <core/mesh/LineElement.hpp>
#include <core/mesh/TriElement.hpp>
#include <core/mesh/TetElement.hpp>

namespace Core {

    PrescribedCurrentDensity::PrescribedCurrentDensity(int element_id, const Eigen::Vector3d& current_density, const std::string& tag)
        : SourceTerm(tag), element_id_(element_id), current_density_(current_density) {}

    void PrescribedCurrentDensity::apply(Eigen::VectorXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const {
        Element* elem = mesh.getElement(element_id_);
        if (!elem) return;

        elem->setOrder(element_order);

        const auto& ref_data = ReferenceElementCache::get(elem->getTypeName(), elem->getNodes().size(), element_order, element_order);
        FEValues fe_values(elem->getGeometry(), element_order, ref_data);

        // **【FIX START】**: Use the comprehensive DOF gathering logic, identical to PhysicsField::getElementDofs
        const auto& vertex_nodes = elem->getNodes();
        const size_t num_vertices = vertex_nodes.size();
        const size_t num_elem_nodes = elem->getNumNodes();

        std::vector<int> base_dofs;
        base_dofs.reserve(num_elem_nodes);
        for (size_t i = 0; i < num_vertices; ++i) {
            base_dofs.push_back(dof_manager.getEquationIndex(vertex_nodes[i]->getId(), var_name));
        }

        if (element_order > 1) {
            if (dynamic_cast<const Core::LineElement*>(elem)) {
                base_dofs.push_back(dof_manager.getEdgeEquationIndex({vertex_nodes[0]->getId(), vertex_nodes[1]->getId()}, var_name));
            } else if (dynamic_cast<const Core::TriElement*>(elem)) {
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto& edge : edges) {
                    base_dofs.push_back(dof_manager.getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, var_name));
                }
            } else if (dynamic_cast<const Core::TetElement*>(elem)) {
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                for (const auto& edge : edges) {
                    base_dofs.push_back(dof_manager.getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, var_name));
                }
            }
        }

        std::vector<int> dofs;
        dofs.reserve(num_elem_nodes * 3); // 3 components for vector potential
        for (int base_dof : base_dofs) {
            if (base_dof != -1) {
                for (int c = 0; c < 3; ++c) dofs.push_back(base_dof + c);
            } else {
                for (int c = 0; c < 3; ++c) dofs.push_back(-1);
            }
        }
        // **【FIX END】**

        Eigen::VectorXd fe_local = Eigen::VectorXd::Zero(dofs.size());

        for (size_t q_p = 0; q_p < fe_values.num_quadrature_points(); ++q_p) {
            fe_values.reinit(q_p);
            const auto& N = fe_values.get_shape_values();
            const double detJ_x_w = fe_values.get_detJ_times_weight();

            for (size_t i = 0; i < num_elem_nodes; ++i) {
                fe_local(i * 3 + 0) += N(i) * current_density_(0) * detJ_x_w;
                fe_local(i * 3 + 1) += N(i) * current_density_(1) * detJ_x_w;
                fe_local(i * 3 + 2) += N(i) * current_density_(2) * detJ_x_w;
            }
        }

        for (size_t i = 0; i < dofs.size(); ++i) {
            if (dofs[i] != -1) {
                F(dofs[i]) += fe_local(i);
            }
        }
    }

} // namespace Core