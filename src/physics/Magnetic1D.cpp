#include "physics/Magnetic1D.hpp"
#include <core/mesh/LineElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"

namespace Physics {

Magnetic1D::Magnetic1D(const Core::Material& material)
    : material_(material) {}

const char* Magnetic1D::getName() const { return "Magnetic Field 1D"; }
const char* Magnetic1D::getVariableName() const { return "MagneticPotential"; }

void Magnetic1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);

    K_.setZero();
    F_.setZero();
    U_.setZero();
}

void Magnetic1D::assemble() {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    F_.setZero();

    const double mu = material_.getProperty("magnetic_permeability");
    const double inv_mu = 1.0 / mu;

    std::vector<Eigen::Triplet<double>> k_triplets;
    auto quadrature_points = Utils::Quadrature::getLineQuadrature(element_order_);

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
        if (line_elem) {
            const auto& vertex_nodes = line_elem->getNodes();
            if (vertex_nodes.size() != 2) continue;

            int order = element_order_;
            size_t num_elem_nodes = order + 1;

            // --- THE FIX IS HERE: Correct DOF Ordering ---
            std::vector<int> dofs(num_elem_nodes);
            if (order == 1) {
                dofs[0] = dof_manager_->getEquationIndex(vertex_nodes[0]->getId(), getVariableName());
                dofs[1] = dof_manager_->getEquationIndex(vertex_nodes[1]->getId(), getVariableName());
            } else if (order == 2) {
                // Order must match shape functions: [Vertex 1, Midpoint, Vertex 2]
                dofs[0] = dof_manager_->getEquationIndex(vertex_nodes[0]->getId(), getVariableName());
                dofs[1] = dof_manager_->getEdgeEquationIndex({vertex_nodes[0]->getId(), vertex_nodes[1]->getId()}, getVariableName());
                dofs[2] = dof_manager_->getEquationIndex(vertex_nodes[1]->getId(), getVariableName());
            } else {
                throw std::runtime_error("Invalid element order.");
            }

            const double h = line_elem->getLength();
            const double detJ = h / 2.0;

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            for(const auto& qp : quadrature_points) {
                Eigen::VectorXd dN_dxi = Utils::ShapeFunctions::getLineShapeFunctionDerivatives(order, qp.point(0));
                Eigen::MatrixXd B = (1.0 / detJ) * dN_dxi.transpose();

                ke_local += B.transpose() * inv_mu * B * qp.weight * detJ;
            }

            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
                     if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                    }
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    logger.info("Assembly for ", getName(), " complete.");
}

}