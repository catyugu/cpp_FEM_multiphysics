#include "physics/Heat3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"
#include <cmath>

namespace Physics {

Heat3D::Heat3D(const Core::Material& material) : material_(material), k_(0.0) {}

const char* Heat3D::getName() const { return "Heat Transfer 3D"; }
const char* Heat3D::getVariableName() const { return "Temperature"; }

void Heat3D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    k_ = material_.getProperty("thermal_conductivity");

    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> Thermal Conductivity (k): ", k_);

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

void Heat3D::assemble(const PhysicsField *coupled_field) {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    M_.setZero();
    applySources(); // Apply source terms to F_

    const Eigen::Matrix3d D_mat = Eigen::Matrix3d::Identity() * k_;
    const double rho_cp = material_.getProperty("density") * material_.getProperty("thermal_capacity");

    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets;

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            // Set the element's mathematical order before creating FEValues
            tet_elem->setOrder(element_order_);

            // 1. Create the FEValues calculator for this element.
            auto fe_values = tet_elem->create_fe_values(element_order_); // Use element_order_ for quad_order too for simplicity

            // 2. Get the correct DOF indices from the centralized function.
            const auto dofs = get_element_dofs(tet_elem);
            const size_t num_elem_nodes = tet_elem->getNumNodes(); // This now correctly reflects the order

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            // 3. Loop over quadrature points using FEValues.
            for(size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);

                // 4. Get pre-calculated values.
                const auto& N = fe_values->get_shape_values();      // Shape function values
                const auto& B = fe_values->get_shape_gradients();   // Shape function gradients in real coordinates (B matrix)
                const double detJ_x_w = fe_values->get_detJ_times_weight(); // det(J) * w_q

                // 5. Compute local matrices.
                // Stiffness matrix: integral( (∇N)^T * D * (∇N) dV )
                ke_local += B.transpose() * D_mat * B * detJ_x_w;
                // Mass matrix: integral( N * rho_cp * N^T dV )
                me_local += N * rho_cp * N.transpose() * detJ_x_w;
            }

            // Assembly into global triplet lists remains the same
            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                        m_triplets.emplace_back(dofs[i], dofs[j], me_local(i, j));
                    }
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());
    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics