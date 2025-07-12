#include "core/sources/VolumetricSource.hpp"
#include "core/DOFManager.hpp"
#include "core/mesh/Mesh.hpp"
#include "core/mesh/Element.hpp"

namespace Core {

    VolumetricSource::VolumetricSource(int element_id, double total_power, const std::string& tag)
        : SourceTerm(tag), element_id_(element_id), total_power_(total_power) {}

    void VolumetricSource::apply(Eigen::MatrixXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name) const {
        const Element* elem = mesh.getElement(element_id_);
        if (!elem) return;

        double nodal_share = total_power_ / elem->getNumNodes();

        for (const auto& node : elem->getNodes()) {
            // The correct var_name is used here, guaranteed by the calling field
            int dof_index = dof_manager.getEquationIndex(node->getId(), var_name);
            if (dof_index != -1) {
                F(dof_index) += nodal_share;
            }
        }
    }

}