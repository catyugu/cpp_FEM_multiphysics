#include "physics/Heat2D.hpp"
#include <core/mesh/TriElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "utils/InterpolationUtilities.hpp"

namespace Physics {
    Heat2D::Heat2D() = default;

    const char *Heat2D::getName() const { return "Heat Transfer 2D"; }
    const char *Heat2D::getVariableName() const { return "Temperature"; }

    void Heat2D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");
    }

    void Heat2D::assemble() {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        M_.setZero();
        applySources();

        std::vector<Eigen::Triplet<double> > k_triplets;
        std::vector<Eigen::Triplet<double> > m_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // 获取材料引用
            const auto& material = getMaterial(elem_ptr);

            auto fe_values = elem_ptr->createFEValues(element_order_);
            fe_values->setAnalysisType(Core::AnalysisType::SCALAR_DIFFUSION);

            const auto dofs = getElementDofs(elem_ptr);
            const auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            // 准备物理场映射以便插值
            std::map<std::string, const Physics::PhysicsField*> physics_fields;
            if (problem_) {
                auto* heat_field = problem_->getField("Temperature");
                auto* voltage_field = problem_->getField("Voltage");
                if (heat_field) physics_fields["Temperature"] = heat_field;
                if (voltage_field) physics_fields["Voltage"] = voltage_field;
            }

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));

                const auto &N = fe_values->get_shape_values();
                const auto &B = fe_values->getBMatrix();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                // ========= 新增：逐积分点插值计算材料属性 =========
                std::vector<std::string> variable_names = {"Temperature"};
                auto interpolated_vars = Utils::InterpolationUtilities::interpolateAtQuadraturePoint(
                    elem_ptr, N, variable_names, physics_fields);

                const double k = material.getPropertyAtQuadraturePoint("thermal_conductivity", interpolated_vars);
                const double rho = material.getPropertyAtQuadraturePoint("density", interpolated_vars);
                const double cp = material.getPropertyAtQuadraturePoint("thermal_capacity", interpolated_vars);
                const double rho_cp = rho * cp;
                const Eigen::Matrix2d D_mat = Eigen::Matrix2d::Identity() * k;
                // ================================================

                ke_local += B.transpose() * D_mat * B * detJ_x_w;
                me_local += N * rho_cp * N.transpose() * detJ_x_w;
            }

            for (size_t i = 0; i < static_cast<size_t>(num_elem_nodes); ++i) {
                for (size_t j = 0; j < static_cast<size_t>(num_elem_nodes); ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
                        m_triplets.emplace_back(dofs[i], dofs[j], me_local(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
                    }
                }
            }
        }
        K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
        M_.setFromTriplets(m_triplets.begin(), m_triplets.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
}
