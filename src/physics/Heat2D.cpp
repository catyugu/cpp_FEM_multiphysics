#include "physics/Heat2D.hpp"
#include <core/mesh/TriElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"

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

    void Heat2D::assemble(const PhysicsField *coupled_field) {
        auto &logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        M_.setZero();
        applySources();

        std::vector<Eigen::Triplet<double> > k_triplets;
        std::vector<Eigen::Triplet<double> > m_triplets;

        for (const auto &elem_ptr: mesh_->getElements()) {
            elem_ptr->setOrder(element_order_);
            
            // --- NEW: Get material for the current element ---
            const auto& material = getMaterial(elem_ptr);
            const double k = material.getProperty("thermal_conductivity");
            const double rho_cp = material.getProperty("density") * material.getProperty("thermal_capacity");
            const Eigen::Matrix2d D_mat = Eigen::Matrix2d::Identity() * k;
            // ------------------------------------------------

            auto fe_values = elem_ptr->createFEValues(element_order_);

            // 新增：设置分析类型为标量扩散问题，自动构建B矩阵
            fe_values->setAnalysisType(Core::AnalysisType::SCALAR_DIFFUSION);

            const auto dofs = getElementDofs(elem_ptr);
            const auto num_elem_nodes = static_cast<Eigen::Index>(elem_ptr->getNumNodes());

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(static_cast<int>(q_p));

                const auto &N = fe_values->get_shape_values();
                // 直接获取预构建的B矩阵（梯度矩阵）
                const auto &B = fe_values->getBMatrix();
                const double detJ_x_w = fe_values->get_detJ_times_weight();

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
