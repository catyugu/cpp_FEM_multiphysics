#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/LineElement.hpp>
#include "core/FEValues.hpp" // Use the FEValues calculator
#include "core/sources/SourceTerm.hpp"

namespace Physics {

Heat1D::Heat1D(const Core::Material& material) : material_(material) {}

const char* Heat1D::getName() const { return "Heat Transfer 1D"; }
const char* Heat1D::getVariableName() const { return "Temperature"; }

void Heat1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);

    F_.setZero();
    U_.setZero();
    U_prev_.setZero();
}

void Heat1D::assemble(const PhysicsField *coupled_field) {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    M_.setZero();
    applySources();

    const double k = material_.getProperty("thermal_conductivity");
    const double rho_cp = material_.getProperty("density") * material_.getProperty("thermal_capacity");

    std::vector<Eigen::Triplet<double>> k_triplets, m_triplets;

    for (const auto& elem_ptr : mesh_->getElements()) {
        if (auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr)) {
            line_elem->setOrder(element_order_);

            auto fe_values = line_elem->create_fe_values(element_order_);
            const auto dofs = get_element_dofs(line_elem);
            const size_t num_elem_nodes = line_elem->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for(size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);
                const auto& N = fe_values->get_shape_values();
                const auto& B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                ke_local += B.transpose() * k * B * detJ_x_w;
                me_local += N * rho_cp * N.transpose() * detJ_x_w;
            }

            for(size_t r = 0; r < num_elem_nodes; ++r) {
                for(size_t c = 0; c < num_elem_nodes; ++c) {
                    if (dofs[r] != -1 && dofs[c] != -1) {
                        k_triplets.emplace_back(dofs[r], dofs[c], ke_local(r,c));
                        m_triplets.emplace_back(dofs[r], dofs[c], me_local(r,c));
                    }
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());
    logger.info("Assembly for ", getName(), " complete.");
}

}