#include "physics/Magnetic1D.hpp"
#include <core/mesh/LineElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp" // Use the FEValues calculator

namespace Physics {

Magnetic1D::Magnetic1D(const Core::Material& material)
    : material_(material) {}

const char* Magnetic1D::getName() const { return "Magnetic Field 1D"; }
const char* Magnetic1D::getVariableName() const { return "MagneticPotential"; }

void Magnetic1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);

    K_.setZero();
    F_.setZero();
    U_.setZero();
}

void Magnetic1D::assemble(const PhysicsField *coupled_field) {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    F_.setZero();

    const double inv_mu = 1.0 / material_.getProperty("magnetic_permeability");

    std::vector<Eigen::Triplet<double>> k_triplets;

    for (const auto& elem_ptr : mesh_->getElements()) {
        if (auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr)) {
            line_elem->setOrder(element_order_);

            auto fe_values = line_elem->create_fe_values(element_order_);
            const auto dofs = getElementDofs(line_elem);
            const size_t num_elem_nodes = line_elem->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            for(size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);
                const auto& B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                ke_local += B.transpose() * inv_mu * B * detJ_x_w;
            }

            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
                     if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                    }
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    logger.info("Assembly for ", getName(), " complete.");
}

}