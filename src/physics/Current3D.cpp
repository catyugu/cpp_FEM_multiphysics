#include "physics/Current3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "core/ReferenceElement.hpp"
#include "utils/Exceptions.hpp"
#include <cmath>
#include "physics/Heat3D.hpp"

namespace Physics {
    Current3D::Current3D() {
    }

    const char *Current3D::getName() const { return "Current 3D"; }
    const char *Current3D::getVariableName() const { return "Voltage"; }

    void Current3D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");

        size_t num_eq = dof_manager_->getNumEquations();
        K_.resize(num_eq, num_eq);
        F_.resize(num_eq, 1);
        U_.resize(num_eq, 1);
        U_prev_.resize(num_eq, 1);
        F_.setZero();
        U_.setZero();
        U_prev_.setZero();
    }

    // src/physics/Current3D.cpp

    void Current3D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();

        const Physics::Heat3D *heat_field = nullptr;
        if (coupled_field) {
            heat_field = dynamic_cast<const Physics::Heat3D *>(coupled_field);
        }
        const Eigen::VectorXd &heat_solution = heat_field ? heat_field->getSolution() : U_;

        std::vector<Eigen::Triplet<double> > triplet_list;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // --- NEW: Get material for the current element ---
            const auto& material = getMaterial(elem_ptr);
            // ------------------------------------------------

            auto fe_values = elem_ptr->createFEValues(element_order_);

            const auto dofs = getElementDofs(elem_ptr);
            const size_t num_elem_nodes = elem_ptr->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            const auto heat_dofs = heat_field ? heat_field->getElementDofs(elem_ptr) : std::vector<int>();
            Eigen::VectorXd nodal_temperatures(heat_dofs.size());
            if (heat_field) {
                for (size_t k = 0; k < heat_dofs.size(); ++k) {
                    nodal_temperatures(k) = (heat_dofs[k] != -1) ? heat_solution(heat_dofs[k]) : 293.15;
                }
            }

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);
                const auto &N = fe_values->get_shape_values();
                const auto &B = fe_values->get_shape_gradients();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                double sigma;
                if (heat_field) {
                    double temp_at_qp = N.transpose() * nodal_temperatures;
                    sigma = material.getProperty("electrical_conductivity", {{"Temperature", temp_at_qp}});
                } else {
                    sigma = material.getProperty("electrical_conductivity");
                }

                const Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * sigma;
                ke_local += B.transpose() * D * B * detJ_x_w;
            }

            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
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
