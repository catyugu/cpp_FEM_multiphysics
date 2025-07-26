#include "physics/Current2D.hpp"
#include <core/mesh/TriElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp" // Use the FEValues calculator

namespace Physics {

    Current2D::Current2D(const Core::Material& material)
        : material_(material) {}

    const char* Current2D::getName() const { return "Current 2D"; }
    const char* Current2D::getVariableName() const { return "Voltage"; }

    void Current2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
        mesh_ = &mesh;
        dof_manager_ = &dof_manager;
        auto& logger = Utils::Logger::instance();
        logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

        size_t num_eq = dof_manager_->getNumEquations();

        K_.resize(num_eq, num_eq);
        F_.resize(num_eq, 1);
        U_.resize(num_eq, 1);
        F_.setZero();
        U_.setZero();
    }


    void Current2D::assemble(const PhysicsField *coupled_field) {
        auto& logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero();

        const double sigma = material_.getProperty("electrical_conductivity");
        const Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * sigma;

        std::vector<Eigen::Triplet<double>> triplet_list;

        for (const auto& elem_ptr : mesh_->getElements()) {
            if (auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr)) {
                // Set the element's mathematical order before creating FEValues
                tri_elem->setOrder(element_order_);

                auto fe_values = tri_elem->create_fe_values(element_order_);
                const auto dofs = get_element_dofs(tri_elem);
                const size_t num_elem_nodes = tri_elem->getNumNodes();

                Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

                for(size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                    fe_values->reinit(q_p);
                    const auto& B = fe_values->get_shape_gradients();
                    const double detJ_x_w = fe_values->get_detJ_times_weight();

                    ke_local += B.transpose() * D * B * detJ_x_w;
                }

                for (size_t i=0; i < num_elem_nodes; ++i) {
                    for (size_t j=0; j < num_elem_nodes; ++j) {
                        if (dofs[i] != -1 && dofs[j] != -1) {
                            triplet_list.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                        }
                    }
                }
            }
        }
        K_.setFromTriplets(triplet_list.begin(), triplet_list.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
}