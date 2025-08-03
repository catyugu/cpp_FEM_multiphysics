// src/post/MagneticFieldCalculator.cpp

#include "post/MagneticFieldCalculator.hpp"
#include "core/Problem.hpp"
#include "physics/PhysicsField.hpp"
#include "core/mesh/Element.hpp"
#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"
#include "utils/SimpleLogger.hpp"

namespace Post {

    const char* MagneticFieldCalculator::getName() const {
        return "MagneticField_B";
    }

    PostProcessingResult MagneticFieldCalculator::compute_derived_quantities(const Core::Problem& problem) const {
        auto& logger = Utils::Logger::instance();
        logger.info("Post-processing: Calculating Magnetic Field (B)...");

        // 1. 获取所需的物理场和数据
        const auto* magnetic_field = problem.getField("MagneticVectorPotential");
        if (!magnetic_field) {
            throw Exception::ConfigurationException("MagneticFieldCalculator requires a 'MagneticVectorPotential' field.");
        }

        const auto& mesh = problem.getMesh();
        const auto& A_solution = magnetic_field->getSolution();

        PostProcessingResult result;
        result.name = getName();
        result.dimension = 3; // B 是一个三维矢量
        result.data.resize(mesh.getElements().size());

        // 2. 遍历网格中的所有单元
        for (size_t i = 0; i < mesh.getElements().size(); ++i) {
            Core::Element* elem = mesh.getElements()[i];
            elem->setOrder(magnetic_field->getElementOrder());

            const auto& ref_data = Core::ReferenceElementCache::get(elem->getTypeName(), elem->getNodes().size(), magnetic_field->getElementOrder(), magnetic_field->getElementOrder());
            Core::FEValues fe_values(elem->getGeometry(), magnetic_field->getElementOrder(), ref_data);

            const auto dofs = magnetic_field->getElementDofs(elem);
            const size_t num_elem_nodes = elem->getNumNodes();
            const int num_components = magnetic_field->getNumComponents();

            // 提取该单元的节点上的磁矢量势 A 的值 (Ax, Ay, Az, ...)
            Eigen::VectorXd nodal_A_values(dofs.size());
            for (size_t j = 0; j < dofs.size(); ++j) {
                nodal_A_values(j) = (dofs[j] != -1) ? A_solution(dofs[j]) : 0.0;
            }

            result.data[i].resize(fe_values.num_quadrature_points());
            // 3. 遍历每个单元的所有积分点 (高斯点)
            for (size_t q_p = 0; q_p < fe_values.num_quadrature_points(); ++q_p) {
                fe_values.reinit(q_p);
                const auto& grad_N = fe_values.get_shape_gradients(); // 获取形函数梯度 ∇N

                // 构建旋度算子矩阵 B_curl
                Eigen::MatrixXd B_curl(3, num_elem_nodes * num_components);
                B_curl.setZero();
                for(size_t node_idx = 0; node_idx < num_elem_nodes; ++node_idx) {
                    double dN_dx = grad_N(0, node_idx);
                    double dN_dy = grad_N(1, node_idx);
                    double dN_dz = grad_N(2, node_idx);

                    // 填充 B_curl 矩阵，这与 Magnetic3D::assemble 中的逻辑完全一致
                    B_curl(0, node_idx*3 + 1) = -dN_dz; B_curl(0, node_idx*3 + 2) =  dN_dy;
                    B_curl(1, node_idx*3 + 0) =  dN_dz; B_curl(1, node_idx*3 + 2) = -dN_dx;
                    B_curl(2, node_idx*3 + 0) = -dN_dy; B_curl(2, node_idx*3 + 1) =  dN_dx;
                }

                // 计算积分点处的磁场 B = B_curl * nodal_A
                result.data[i][q_p] = B_curl * nodal_A_values;
            }
        }

        logger.info("Magnetic Field (B) calculation complete.");
        return result;
    }

} // namespace Post