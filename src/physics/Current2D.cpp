#include "physics/Current2D.hpp"
#include <core/mesh/TriElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/SimpleLogger.hpp"
#include <vector>
#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"
#include <cmath> // For std::abs

namespace Physics {

    Current2D::Current2D(const Core::Material& material)
        : material_(material) {}

    const char* Current2D::getName() const { return "Electromagnetics 2D"; }
    const char* Current2D::getVariableName() const { return "Voltage"; }

    void Current2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
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


    void Current2D::assemble() {
        auto& logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero();

        std::vector<Eigen::Triplet<double>> triplet_list;
        auto quadrature_points = Utils::Quadrature::getTriangleQuadrature(element_order_);

        for (const auto& elem_ptr : mesh_->getElements()) {
            auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
            if (tri_elem) {
                const auto& vertex_nodes = tri_elem->getNodes();
                if (vertex_nodes.size() != 3) continue;

                int order = element_order_;
                size_t num_elem_nodes = (order + 1) * (order + 2) / 2;

                Eigen::Matrix<double, 2, Eigen::Dynamic> all_node_coords(2, num_elem_nodes);
                std::vector<int> dofs(num_elem_nodes);

                for(size_t i = 0; i < 3; ++i) {
                    const auto& coords = vertex_nodes[i]->getCoords();
                    all_node_coords(0, i) = coords[0];
                    all_node_coords(1, i) = coords[1];
                    dofs[i] = dof_manager_->getEquationIndex(vertex_nodes[i]->getId(), getVariableName());
                }

                if (order > 1) {
                    int edge_node_idx = 3;
                    const std::vector<std::pair<int, int>> edges = {{0,1}, {1,2}, {2,0}};
                    for (const auto& edge : edges) {
                        all_node_coords.col(edge_node_idx) = (all_node_coords.col(edge.first) + all_node_coords.col(edge.second)) * 0.5;
                        dofs[edge_node_idx] = dof_manager_->getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName());
                        edge_node_idx++;
                    }
                }

                double local_sigma = material_.getProperty("electrical_conductivity");
                Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);
                Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * local_sigma;

                for(const auto& qp : quadrature_points) {
                    Eigen::MatrixXd dN_dnat = Utils::ShapeFunctions::getTriShapeFunctionDerivatives(order, qp.point(0), qp.point(1));
                    Eigen::Matrix2d J = all_node_coords * dN_dnat;
                    double detJ = std::abs(J.determinant());
                    Eigen::MatrixXd B = (dN_dnat * J.inverse()).transpose();

                    ke_local += B.transpose() * D * B * qp.weight * detJ;
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