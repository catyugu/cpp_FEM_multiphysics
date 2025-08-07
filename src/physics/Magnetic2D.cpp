#include "physics/Magnetic2D.hpp"
#include <core/mesh/TriElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"

namespace Physics {
    Magnetic2D::Magnetic2D() = default;

    const char *Magnetic2D::getName() const { return "Magnetic Field 2D"; }
    const char *Magnetic2D::getVariableName() const { return "MagneticPotential"; }

    void Magnetic2D::setup(Core::Problem& problem, Core::Mesh &mesh, Core::DOFManager &dof_manager) {
        // Call the base class setup
        PhysicsField::setup(problem, mesh, dof_manager);
        
        auto &logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh.");
    }

    void Magnetic2D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero();

        std::vector<Eigen::Triplet<double>> k_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // --- NEW: Get material for the current element ---
            const auto& material = getMaterial(elem_ptr);
            const double inv_mu = 1.0 / material.getProperty("magnetic_permeability");
            // ------------------------------------------------

            auto fe_values = elem_ptr->createFEValues(element_order_);

            // 新增：设置分析类型为标量扩散问题（2D磁势问题类似扩散）
            fe_values->setAnalysisType(Core::AnalysisType::SCALAR_DIFFUSION);

            const auto dofs = getElementDofs(elem_ptr);
            auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for (Eigen::Index q_p = 0; q_p < static_cast<Eigen::Index>(fe_values->num_quadrature_points()); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));

                // 直接获取预构建的B矩阵（梯度矩阵）
                const auto &B = fe_values->getBMatrix();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

                // 使用2x2的单位矩阵，而不是标量
                const Eigen::Matrix2d D_mat = Eigen::Matrix2d::Identity() * inv_mu;
                ke_local += B.transpose() * D_mat * B * detJ_x_w;
            }

            for (Eigen::Index i = 0; i < num_elem_nodes; ++i) {
                for (Eigen::Index j = 0; j < num_elem_nodes; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(static_cast<int>(dofs[i]), static_cast<int>(dofs[j]), ke_local(i, j));
                    }
                }
            }
        }
        K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
} // namespace Physics
