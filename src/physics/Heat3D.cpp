#include "physics/Heat3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"

namespace Physics {
    Heat3D::Heat3D() : k_(0.0) {
    }

    const char *Heat3D::getName() const { return "Heat Transfer 3D"; }
    const char *Heat3D::getVariableName() const { return "Temperature"; }

    void Heat3D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        // Note: We can no longer set a single 'k_' here, as it can vary per element.
        // We will fetch it inside the assembly loop instead.
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");
    }

    void Heat3D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        M_.setZero();
        applySources();

        std::vector<Eigen::Triplet<double> > k_triplets;
        std::vector<Eigen::Triplet<double> > m_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // --- NEW: Get material for the current element ---
            const auto& material = getMaterial(elem_ptr);
            const double k = material.getProperty("thermal_conductivity");
            const double rho_cp = material.getProperty("density") * material.getProperty("thermal_capacity");
            const Eigen::Matrix3d D_mat = Eigen::Matrix3d::Identity() * k;
            // ------------------------------------------------

            auto fe_values = elem_ptr->createFEValues(element_order_);

            const auto dofs = getElementDofs(elem_ptr);
            const auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for (int q_p = 0; q_p < static_cast<int>(fe_values->num_quadrature_points()); ++q_p) {
                fe_values->reinit(q_p);

                const auto &N = fe_values->get_shape_values();
                const auto &B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                ke_local += B.transpose() * D_mat * B * detJ_x_w;
                me_local += N * rho_cp * N.transpose() * detJ_x_w;
            }

            for (Eigen::Index i = 0; i < num_elem_nodes; ++i) {
                for (Eigen::Index j = 0; j < num_elem_nodes; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(static_cast<int>(dofs[i]), static_cast<int>(dofs[j]), ke_local(i, j));
                        m_triplets.emplace_back(static_cast<int>(dofs[i]), static_cast<int>(dofs[j]), me_local(i, j));
                    }
                }
            }
        }
        K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
        M_.setFromTriplets(m_triplets.begin(), m_triplets.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
} // namespace Physics
