#include "physics/Magnetic3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"

namespace Physics {
    Magnetic3D::Magnetic3D() = default;

    const char *Magnetic3D::getName() const { return "Magnetic Field 3D"; }
    const char *Magnetic3D::getVariableName() const { return "MagneticVectorPotential"; }

    void Magnetic3D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");
    }

    void Magnetic3D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero(); // Sources are applied via addSource, not here

        std::vector<Eigen::Triplet<double> > triplet_list;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // --- NEW: Get material for the current element ---
            const auto& material = getMaterial(elem_ptr);
            // 修复：正确使用磁导率而不是其倒数
            const double mu = material.getProperty("magnetic_permeability");
            const double inv_mu = 1.0 / mu;
            // ------------------------------------------------

            auto fe_values = elem_ptr->createFEValues(element_order_);

            const auto dofs = getElementDofs(elem_ptr);
            auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());
            auto num_components = static_cast<Eigen::Index>(getNumComponents());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes * num_components,
                                                             num_elem_nodes * num_components);

            for (Eigen::Index q_p = 0; q_p < static_cast<Eigen::Index>(fe_values->num_quadrature_points()); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));
                const auto &grad_N = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                Eigen::MatrixXd B_curl(3, num_elem_nodes * 3);
                B_curl.setZero();
                for (Eigen::Index i = 0; i < num_elem_nodes; ++i) {
                    double dN_dx = grad_N(0, i);
                    double dN_dy = grad_N(1, i);
                    double dN_dz = grad_N(2, i);

                    B_curl(0, i * 3 + 1) = -dN_dz;
                    B_curl(0, i * 3 + 2) = dN_dy;
                    B_curl(1, i * 3 + 0) = dN_dz;
                    B_curl(1, i * 3 + 2) = -dN_dx;
                    B_curl(2, i * 3 + 0) = -dN_dy;
                    B_curl(2, i * 3 + 1) = dN_dx;
                }
            ke_local += B_curl.transpose() * inv_mu * B_curl * detJ_x_w;
        }

            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(dofs.size()); ++i) {
                for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(dofs.size()); ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        triplet_list.emplace_back(static_cast<int>(dofs[i]), static_cast<int>(dofs[j]), ke_local(i, j));
                    }
                }
            }
        }
        K_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
} // namespace Physics
