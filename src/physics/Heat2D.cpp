#include "physics/Heat2D.hpp"
#include <core/mesh/TriElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp" // Use the FEValues calculator
#include "core/sources/SourceTerm.hpp"

namespace Physics {

Heat2D::Heat2D(const Core::Material& material) : material_(material) {}

const char* Heat2D::getName() const { return "Heat Transfer 2D"; }
const char* Heat2D::getVariableName() const { return "Temperature"; }

void Heat2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq,1); F_.setZero();
    U_.resize(num_eq,1); U_.setZero();
    U_prev_.resize(num_eq,1); U_prev_.setZero();
}

void Heat2D::assemble(const PhysicsField *coupled_field) {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    M_.setZero();
    applySources(); // Apply source terms to F_

    const double k_therm = material_.getProperty("thermal_conductivity");
    const double rho_cp = material_.getProperty("density") * material_.getProperty("thermal_capacity");
    const Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * k_therm;

    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets;

    for (const auto& elem_ptr : mesh_->getElements()) {
        if (auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr)) {
            // Set the element's mathematical order before creating FEValues
            tri_elem->setOrder(element_order_);

            // 1. Create the FEValues calculator for this element.
            auto fe_values = tri_elem->create_fe_values(element_order_);

            // 2. Get the correct DOF indices from the centralized function.
            const auto dofs = getElementDofs(tri_elem);
            const size_t num_elem_nodes = tri_elem->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            // 3. Loop over quadrature points.
            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);

                // 4. Get pre-calculated values.
                const auto& N = fe_values->get_shape_values();
                const auto& B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                // 5. Compute local matrices.
                ke_local += B.transpose() * D * B * detJ_x_w;
                me_local += N * (rho_cp * N.transpose()) * detJ_x_w;
            }

            // 6. Add to global triplets.
            for (size_t r = 0; r < num_elem_nodes; ++r) {
                for (size_t c = 0; c < num_elem_nodes; ++c) {
                    if (dofs[r] != -1 && dofs[c] != -1) {
                        k_triplets.emplace_back(dofs[r], dofs[c], ke_local(r, c));
                        m_triplets.emplace_back(dofs[r], dofs[c], me_local(r, c));
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