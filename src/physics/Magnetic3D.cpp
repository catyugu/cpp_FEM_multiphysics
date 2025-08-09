#include "physics/Magnetic3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "utils/InterpolationUtilities.hpp"

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

    void Magnetic3D::assemble() {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero(); // Sources are applied via addSource, not here

        std::vector<Eigen::Triplet<double> > triplet_list;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // 获取材料引用
            const auto& material = getMaterial(elem_ptr);

            auto fe_values = elem_ptr->createFEValues(element_order_);
            fe_values->setAnalysisType(Core::AnalysisType::VECTOR_CURL);

            const auto dofs = getElementDofs(elem_ptr);
            auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());
            auto num_components = static_cast<Eigen::Index>(getNumComponents());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes * num_components,
                                                             num_elem_nodes * num_components);

            // 准备物理场映射以便插值
            std::map<std::string, const Physics::PhysicsField*> physics_fields;
            if (problem_) {
                // 获取可能影响磁导率的场（如温度场、磁场强度等）
                auto* heat_field = problem_->getField("Temperature");
                auto* magnetic_field = problem_->getField("MagneticVectorPotential");
                if (heat_field) physics_fields["Temperature"] = heat_field;
                if (magnetic_field) physics_fields["MagneticVectorPotential"] = magnetic_field;
            }

            for (Eigen::Index q_p = 0; q_p < static_cast<Eigen::Index>(fe_values->num_quadrature_points()); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));

                const auto &B_curl = fe_values->getBMatrix();
                const auto &N = fe_values->get_shape_values();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                // ========= 新增：逐积分点插值计算材料属性 =========
                // 在当前积分点插��变量值
                std::vector<std::string> variable_names = {"Temperature", "MagneticVectorPotential"};
                auto interpolated_vars = Utils::InterpolationUtilities::interpolateAtQuadraturePoint(
                    elem_ptr, N, variable_names, physics_fields);

                // 使用插值得到的变量值计算材料属性
                const double mu = material.getPropertyAtQuadraturePoint("magnetic_permeability", interpolated_vars);
                const double inv_mu = 1.0 / mu;
                // ================================================

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
