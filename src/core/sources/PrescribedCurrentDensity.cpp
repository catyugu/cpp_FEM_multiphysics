#include "core/sources/PrescribedCurrentDensity.hpp"
#include "core/mesh/Element.hpp"
#include "core/FEValues.hpp"
#include "physics/PhysicsField.hpp"

namespace Core {

    PrescribedCurrentDensity::PrescribedCurrentDensity(int element_id, const Eigen::Vector3d& current_density, const std::string& tag)
        : SourceTerm(tag), element_id_(element_id), current_density_(current_density) {}

    void PrescribedCurrentDensity::apply(Eigen::VectorXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const {
        Element* elem = mesh.getElement(element_id_);
        if (!elem) return;

        elem->setOrder(element_order);

        const auto& ref_data = ReferenceElementCache::get(elem->getTypeName(), elem->getNodes().size(), element_order, element_order);
        FEValues fe_values(elem->getGeometry(), element_order, ref_data);
        // This is a stand-in for a proper get_element_dofs from a PhysicsField context
        // It manually reconstructs the DOF list for all components
        std::vector<int> dofs;
        const size_t num_elem_nodes = elem->getNumNodes();
        dofs.reserve(num_elem_nodes * 3);

        std::vector<int> base_dofs;
        base_dofs.reserve(num_elem_nodes);
        for(const auto& node : elem->getNodes()) {
            base_dofs.push_back(dof_manager.getEquationIndex(node->getId(), var_name));
        }

        for (int base_dof : base_dofs) {
            if (base_dof != -1) {
                for (int c = 0; c < 3; ++c) {
                    dofs.push_back(base_dof + c);
                }
            } else {
                for (int c = 0; c < 3; ++c) {
                    dofs.push_back(-1);
                }
            }
        }

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