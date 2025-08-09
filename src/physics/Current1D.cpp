#include "physics/Current1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/LineElement.hpp>
#include "core/FEValues.hpp"
#include "core/ReferenceElement.hpp"
#include "utils/InterpolationUtilities.hpp"

namespace Physics {
    Current1D::Current1D() {
    }

    const char *Current1D::getName() const { return "Current 1D"; }
    const char *Current1D::getVariableName() const { return "Voltage"; }

    void Current1D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");
    }

    void Current1D::assemble() {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero();

        std::vector<Eigen::Triplet<double>> k_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // 获取材料引用
            const auto& material = getMaterial(elem_ptr);

            auto fe_values = elem_ptr->createFEValues(element_order_);
            fe_values->setAnalysisType(Core::AnalysisType::SCALAR_DIFFUSION);

            const auto dofs = getElementDofs(elem_ptr);
            const size_t num_elem_nodes = elem_ptr->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            // 准备物理场映射以便插值
            std::map<std::string, const Physics::PhysicsField*> physics_fields;
            if (problem_) {
                auto* heat_field = problem_->getField("Temperature");
                auto* voltage_field = problem_->getField("Voltage");
                if (heat_field) physics_fields["Temperature"] = heat_field;
                if (voltage_field) physics_fields["Voltage"] = voltage_field;
            }

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);
                const auto &B = fe_values->getBMatrix();
                const auto &N = fe_values->get_shape_values();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                // ========= 新增：逐积分点插值计算材料属性 =========
                std::vector<std::string> variable_names = {"Temperature"};
                auto interpolated_vars = Utils::InterpolationUtilities::interpolateAtQuadraturePoint(
                    elem_ptr, N, variable_names, physics_fields);

                const double local_sigma = material.getPropertyAtQuadraturePoint("electrical_conductivity", interpolated_vars);
                // ================================================

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
