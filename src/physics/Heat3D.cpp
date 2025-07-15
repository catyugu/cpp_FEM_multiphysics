#include "physics/Heat3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"
#include "utils/Exceptions.hpp"
#include "utils/ShapeFunctions.hpp"
#include <cmath> // Required for std::abs

namespace Physics {

Heat3D::Heat3D(const Core::Material& material) : material_(material), k_(0.0) {}

const char* Heat3D::getName() const { return "Heat Transfer 3D"; }
const char* Heat3D::getVariableName() const { return "Temperature"; }

void Heat3D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    k_ = material_.getProperty("thermal_conductivity");

    auto& logger = Utils::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> Thermal Conductivity (k): ", k_);

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

void Heat3D::assemble() {
    auto& logger = Utils::Logger::instance();
    logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

    K_.setZero();
    M_.setZero();

    const Eigen::Matrix3d D_mat = Eigen::Matrix3d::Identity() * k_;
    const double rho_cp = material_.getProperty("density") * material_.getProperty("specific_heat");

    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets;
    auto quadrature_points = Utils::Quadrature::getTetrahedronQuadrature(element_order_);

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            const auto& vertex_nodes = tet_elem->getNodes();
            if (vertex_nodes.size() != 4) continue;

            int order = element_order_;
            // --- FIX IS HERE: Direct calculation for the number of nodes ---
            size_t num_elem_nodes = (order + 1) * (order + 2) * (order + 3) / 6;

            Eigen::Matrix<double, 3, Eigen::Dynamic> all_node_coords(3, num_elem_nodes);
            std::vector<int> dofs(num_elem_nodes);

            for(size_t i = 0; i < 4; ++i) {
                const auto& coords = vertex_nodes[i]->getCoords();
                all_node_coords(0, i) = coords[0];
                all_node_coords(1, i) = coords[1];
                all_node_coords(2, i) = coords[2];
                dofs[i] = dof_manager_->getEquationIndex(vertex_nodes[i]->getId(), getVariableName());
            }

            if (order > 1) {
                int edge_node_idx = 4;
                const std::vector<std::pair<int, int>> edges = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
                for (const auto& edge : edges) {
                    all_node_coords.col(edge_node_idx) = (all_node_coords.col(edge.first) + all_node_coords.col(edge.second)) * 0.5;
                    dofs[edge_node_idx] = dof_manager_->getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName());
                    edge_node_idx++;
                }
            }

            Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
            Eigen::MatrixXd me_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

            for(const auto& qp : quadrature_points) {
                Eigen::VectorXd N = Utils::ShapeFunctions::getTetShapeFunctions(order, qp.point(0), qp.point(1), qp.point(2));
                Eigen::MatrixXd dN_dnat = Utils::ShapeFunctions::getTetShapeFunctionDerivatives(order, qp.point(0), qp.point(1), qp.point(2));

                Eigen::Matrix3d J = all_node_coords * dN_dnat;
                double detJ = std::abs(J.determinant());

                Eigen::MatrixXd B = (dN_dnat * J.inverse()).transpose();

                ke_local += B.transpose() * D_mat * B * qp.weight * detJ;
                me_local += N * rho_cp * N.transpose() * qp.weight * detJ;
            }

            for (size_t i = 0; i < num_elem_nodes; ++i) {
                for (size_t j = 0; j < num_elem_nodes; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                        m_triplets.emplace_back(dofs[i], dofs[j], me_local(i, j));
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