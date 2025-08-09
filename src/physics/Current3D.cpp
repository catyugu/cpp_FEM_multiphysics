#include "physics/Current3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "utils/InterpolationUtilities.hpp"
#include "physics/Heat3D.hpp"

namespace Physics {
    Current3D::Current3D() = default;

    const char *Current3D::getName() const { return "Current 3D"; }
    const char *Current3D::getVariableName() const { return "Voltage"; }

    void Current3D::setup(Core::Problem &problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);

        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");

        const auto num_eq = static_cast<Eigen::Index>(dof_manager_->getNumEquations());
        K_.resize(num_eq, num_eq);
        F_.resize(num_eq, 1);
        U_.resize(num_eq, 1);
        U_prev_.resize(num_eq, 1);
        F_.setZero();
        U_.setZero();
        U_prev_.setZero();
    }

    // src/physics/Current3D.cpp

    void Current3D::assemble() {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();

        std::vector<Eigen::Triplet<double> > triplet_list;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);

            // 获取材料引用
            const auto &material = getMaterial(elem_ptr);

            auto fe_values = elem_ptr->createFEValues(element_order_);
            fe_values->setAnalysisType(Core::AnalysisType::SCALAR_DIFFUSION);

            const auto dofs = getElementDofs(elem_ptr);
            const auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            // 准备物理场映射以便插值
            std::map<std::string, const Physics::PhysicsField*> physics_fields;
            if (problem_) {
                // 获取温度场和电压场用于插值
                auto* heat_field = problem_->getField("Temperature");
                auto* voltage_field = problem_->getField("Voltage");
                if (heat_field) physics_fields["Temperature"] = heat_field;
                if (voltage_field) physics_fields["Voltage"] = voltage_field;
            }

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));
                const auto &B = fe_values->getBMatrix();
                const auto &N = fe_values->get_shape_values();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                // ========= 新增：逐积分点插值计算材料属性 =========
                // 在当前积分点插值变量值
                std::vector<std::string> variable_names = {"Temperature"};
                auto interpolated_vars = Utils::InterpolationUtilities::interpolateAtQuadraturePoint(
                    elem_ptr, N, variable_names, physics_fields);

                // 使用插值得到的变量值计算材料属性
                const double sigma = material.getPropertyAtQuadraturePoint("electrical_conductivity", interpolated_vars);
                const Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * sigma;
                // ================================================

                ke_local += B.transpose() * D * B * detJ_x_w;
            }

            for (size_t i = 0; i < dofs.size(); ++i) {
                for (size_t j = 0; j < dofs.size(); ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        triplet_list.emplace_back(dofs[i], dofs[j], ke_local(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
                    }
                }
            }
        }
        K_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
} // namespace Physics
