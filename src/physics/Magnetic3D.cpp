#include "physics/Magnetic3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "core/ReferenceElement.hpp"
#include "utils/Exceptions.hpp"
#include <cmath>

namespace Physics {
    Magnetic3D::Magnetic3D(const Core::Material &material) : material_(material) {
    }

    const char *Magnetic3D::getName() const { return "Magnetic Field 3D"; }
    const char *Magnetic3D::getVariableName() const { return "MagneticVectorPotential"; }

    void Magnetic3D::setup(Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        mesh_ = &mesh;
        dof_manager_ = &dof_manager;
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

        size_t num_eq = dof_manager_->getNumEquations();
        K_.resize(num_eq, num_eq);
        F_.resize(num_eq, 1);
        U_.resize(num_eq, 1);
        U_prev_.resize(num_eq, 1);
        F_.setZero();
        U_.setZero();
        U_prev_.setZero();
    }

    void Magnetic3D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();

        const double inv_mu = 1.0 / material_.getProperty("magnetic_permeability");

        std::vector<Eigen::Triplet<double> > triplet_list;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);

            // --- 重构后的代码 ---
            auto fe_values = elem_ptr->createFEValues(element_order_);
            // --------------------

            const auto dofs = getElementDofs(elem_ptr);
            const size_t num_elem_nodes = elem_ptr->getNumNodes();
            const int num_components = getNumComponents();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes * num_components,
                                                             num_elem_nodes * num_components);

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);
                const auto &grad_N = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                Eigen::MatrixXd B_curl(3, num_elem_nodes * 3);
                B_curl.setZero();
                for (size_t i = 0; i < num_elem_nodes; ++i) {
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

            for (size_t i = 0; i < dofs.size(); ++i) {
                for (size_t j = 0; j < dofs.size(); ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        triplet_list.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                    }
                }
            }
        }
        K_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
} // namespace Physics
