#include "physics/Magnetic3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"
#include <cmath>

namespace Physics {

Magnetic3D::Magnetic3D(const Core::Material& material) : material_(material) {}

const char* Magnetic3D::getName() const { return "Magnetic Field 3D"; }
const char* Magnetic3D::getVariableName() const { return "MagneticPotential"; }

void Magnetic3D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);
    F_.setZero();
    U_.setZero();
    U_prev_.setZero();
}

void Magnetic3D::assemble() {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    F_.setZero();

    const double mu = material_.getProperty("magnetic_permeability");
    const double inv_mu = 1.0 / mu;
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * inv_mu;

    std::vector<Eigen::Triplet<double>> triplet_list;

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            // Set the element's mathematical order before creating FEValues
            tet_elem->setOrder(element_order_);

            // 1. Create the FEValues calculator for this element.
            auto fe_values = tet_elem->create_fe_values(element_order_);

            // 2. Get the correct DOF indices from the centralized function.
            const auto dofs = get_element_dofs(tet_elem);
            const size_t num_elem_nodes = tet_elem->getNumNodes();

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            // 3. Loop over quadrature points using FEValues.
            for(size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                fe_values->reinit(q_p);

                // 4. Get pre-calculated values (gradients in REAL coordinates).
                const auto& B = fe_values->get_shape_gradients(); // This is the B matrix (∇N)
                const double detJ_x_w = fe_values->get_detJ_times_weight(); // This is det(J) * w_q

                // 5. Compute local matrices using these values.
                // For Magnetic, it's integral( (∇N)^T * D * (∇N) dV )
                ke_local += B.transpose() * D * B * detJ_x_w;
            }

            // Assembly into global triplet list remains the same
            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
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

} // namespace Physics