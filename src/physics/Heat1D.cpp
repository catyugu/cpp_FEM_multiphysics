#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/LineElement.hpp>
#include "core/FEValues.hpp"
#include "core/sources/SourceTerm.hpp"

namespace Physics {
    Heat1D::Heat1D() = default;

    const char *Heat1D::getName() const { return "Heat Transfer 1D"; }
    const char *Heat1D::getVariableName() const { return "Temperature"; }

    void Heat1D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");
    }

    void Heat1D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        M_.setZero();
        applySources();

        std::vector<Eigen::Triplet<double> > k_triplets, m_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // --- NEW: Get material for the current element ---
            const auto& material = getMaterial(elem_ptr);
            const double k = material.getProperty("thermal_conductivity");
            const double rho_cp = material.getProperty("density") * material.getProperty("thermal_capacity");
            // ------------------------------------------------

            auto fe_values = elem_ptr->createFEValues(element_order_);

            const auto dofs = getElementDofs(elem_ptr);
            const auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));
                const auto &N = fe_values->get_shape_values();
                const auto &B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                ke_local += B.transpose() * k * B * detJ_x_w;
                me_local += N * rho_cp * N.transpose() * detJ_x_w;
            }

            for (size_t r = 0; r < dofs.size(); ++r) {
                for (size_t c = 0; c < dofs.size(); ++c) {
                    if (dofs[r] != -1 && dofs[c] != -1) {
                        k_triplets.emplace_back(dofs[r], dofs[c], ke_local(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)));
                        m_triplets.emplace_back(dofs[r], dofs[c], me_local(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)));
                    }
                }
            }
        }
        K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
        M_.setFromTriplets(m_triplets.begin(), m_triplets.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
} // namespace Physics