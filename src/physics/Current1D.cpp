#include "physics/Current1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/LineElement.hpp>
#include "core/FEValues.hpp"
#include "core/ReferenceElement.hpp"

namespace Physics {
    Current1D::Current1D(const Core::Material &material)
        : material_(material) {
    }

    const char *Current1D::getName() const { return "Current 1D"; }
    const char *Current1D::getVariableName() const { return "Voltage"; }

    void Current1D::setup(Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        mesh_ = &mesh;
        dof_manager_ = &dof_manager;
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
        size_t num_eq = dof_manager_->getNumEquations();
        K_.resize(num_eq, num_eq);
        M_.resize(num_eq, num_eq);
        F_.resize(num_eq, 1);
        U_.resize(num_eq, 1);
        F_.setZero();
        U_.setZero();
    }

    void Current1D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero();

        const double local_sigma = material_.getProperty("electrical_conductivity");
        std::vector<Eigen::Triplet<double>> k_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);

            // --- 重构后的代码 ---
            auto fe_values = elem_ptr->create_fe_values(element_order_);
            // --------------------

            const auto dofs = getElementDofs(elem_ptr);
            const size_t num_elem_nodes = elem_ptr->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);
                const auto &B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                ke_local += B.transpose() * local_sigma * B * detJ_x_w;
            }

            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                    }
                }
            }
        }
        K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
}
