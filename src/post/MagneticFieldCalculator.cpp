#include "post/MagneticFieldCalculator.hpp"
#include "core/Problem.hpp"
#include "physics/PhysicsField.hpp"
#include "core/mesh/Element.hpp"
#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"
#include "utils/SimpleLogger.hpp"

namespace Post {

    const char* MagneticFieldCalculator::getName() const {
        return "MagneticField";
    }

    PostProcessingResult MagneticFieldCalculator::compute_derived_quantities(const Core::Problem& problem) const {
        auto& logger = Utils::Logger::instance();
        logger.info("Post-processing: Calculating Magnetic Field (B)...");

        // 1. 获取所需的物理场和数据
        const auto* magnetic_field = problem.getField("MagneticVectorPotential");
        if (!magnetic_field) {
            throw Exception::ConfigurationException("MagneticFieldCalculator requires a 'MagneticVectorPotential' field.");
        }
        if (magnetic_field->getNumComponents() != 3) {
            throw Exception::ConfigurationException("MagneticVectorPotential field must be a 3-component vector field.");
        }

        const auto& mesh = problem.getMesh();
        const auto& A_solution = magnetic_field->getSolution();

        PostProcessingResult result;
        result.name = getName();
        result.dimension = 3; // B field is a 3D vector
        result.data.resize(mesh.getElements().size());

        // 2. 遍历网格中的所有单元
        for (size_t i = 0; i < mesh.getElements().size(); ++i) {
            Core::Element* elem = mesh.getElements()[i];
            elem->setOrder(magnetic_field->getElementOrder());

            const auto& ref_data = Core::ReferenceElementCache::get(elem->getTypeName(), elem->getNodes().size(), magnetic_field->getElementOrder(), magnetic_field->getElementOrder());
            Core::FEValues fe_values(elem->getGeometry(), magnetic_field->getElementOrder(), ref_data);

            const auto dofs = magnetic_field->getElementDofs(elem);
            const size_t num_elem_nodes = elem->getNumNodes();

            // 3. 提取当前单元的节点解 (A_x, A_y, A_z at each node)
            Eigen::VectorXd nodal_A_values(dofs.size());
            for (size_t j = 0; j < dofs.size(); ++j) {
                nodal_A_values(j) = (dofs[j] != -1) ? A_solution(dofs[j]) : 0.0;
            }

            // 4. 在每个单元的积分点上计算B场
            result.data[i].resize(fe_values.num_quadrature_points());
            for (size_t q_p = 0; q_p < fe_values.num_quadrature_points(); ++q_p) {
                fe_values.reinit(q_p);
                const auto& grad_N = fe_values.get_shape_gradients(); // ∇N in real coordinates

                // 5. Build the B_curl matrix with the CORRECTED signs
                Eigen::MatrixXd B_curl(3, num_elem_nodes * 3);
                B_curl.setZero();
                for(size_t node_idx = 0; node_idx < num_elem_nodes; ++node_idx) {
                    double dN_dx = grad_N(0, node_idx);
                    double dN_dy = grad_N(1, node_idx);
                    double dN_dz = grad_N(2, node_idx);

                    // Correct Bx = dAz/dy - dAy/dz
                    B_curl(0, node_idx*3 + 1) = -dN_dz; B_curl(0, node_idx*3 + 2) =  dN_dy;
                    // Correct By = dAx/dz - dAz/dx
                    B_curl(1, node_idx*3 + 0) =  dN_dz; B_curl(1, node_idx*3 + 2) = -dN_dx;
                    // Correct Bz = dAy/dx - dAx/dy
                    B_curl(2, node_idx*3 + 0) = -dN_dy; B_curl(2, node_idx*3 + 1) =  dN_dx;
                }

                // 6. Calculate B = B_curl * A_nodal_values
                result.data[i][q_p] = B_curl * nodal_A_values;
            }
        }

        logger.info("Magnetic Field calculation complete.");
        return result;
    }

} // namespace Post