#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/LineElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/Quadrature.hpp"
#include "core/sources/SourceTerm.hpp"
#include "utils/ShapeFunctions.hpp"

namespace Physics {

Heat1D::Heat1D(const Core::Material& material) : material_(material) {}

const char* Heat1D::getName() const { return "Heat Transfer 1D"; }
const char* Heat1D::getVariableName() const { return "Temperature"; }

void Heat1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);

    F_.setZero();
    U_.setZero();
    U_prev_.setZero();
}

void Heat1D::assemble() {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    M_.setZero();
    F_.setZero();
    for (const auto& source : source_terms_) {
        source->apply(F_, *dof_manager_, *mesh_, getVariableName());
    }

    const double k = material_.getProperty("thermal_conductivity");
    const double rho = material_.getProperty("density");
    const double cp = material_.getProperty("specific_heat");

    auto quadrature_points = Utils::Quadrature::getLineQuadrature(element_order_);
    std::vector<Eigen::Triplet<double>> k_triplets, m_triplets;

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
                throw std::runtime_error("Unsupported element order.");
            }

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            const double h = line_elem->getLength();
            const double detJ = h / 2.0;

            for(const auto& qp : quadrature_points) {
                Eigen::VectorXd N = Utils::ShapeFunctions::getLineShapeFunctions(order, qp.point(0));
                Eigen::VectorXd dN_dxi = Utils::ShapeFunctions::getLineShapeFunctionDerivatives(order, qp.point(0));

                Eigen::MatrixXd B = (1.0 / detJ) * dN_dxi.transpose();

                ke_local += B.transpose() * k * B * qp.weight * detJ;
                me_local += N * (rho * cp) * N.transpose() * qp.weight * detJ;
            }

            for(size_t r = 0; r < num_elem_nodes; ++r) {
                for(size_t c = 0; c < num_elem_nodes; ++c) {
                    if (dofs[r] != -1 && dofs[c] != -1) {
                        k_triplets.emplace_back(dofs[r], dofs[c], ke_local(r,c));
                        m_triplets.emplace_back(dofs[r], dofs[c], me_local(r,c));
                    }
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());
    logger.info("Assembly for ", getName(), " complete.");
}

}