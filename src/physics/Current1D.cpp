#include "physics/Current1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/LineElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"

namespace Physics {
    Current1D::Current1D(const Core::Material& material)
        : material_(material) {}

    const char* Current1D::getName() const { return "Electromagnetics 1D"; }
    const char* Current1D::getVariableName() const { return "Voltage"; }

    void Current1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
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
  void Current1D::assemble() {
        auto& logger = Utils::Logger::instance();
        logger.info("Assembling system for ", getName(), " using mathematical order ", element_order_);

        K_.setZero();
        F_.setZero();

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
                    throw std::runtime_error("Unsupported element order.");
                }
                // --- End of Fix ---

                double local_sigma = material_.getProperty("electrical_conductivity");
                const double h = line_elem->getLength();
                const double detJ = h / 2.0;

                Eigen::MatrixXd ke_local = Eigen::MatrixXd::Zero(num_elem_nodes, num_elem_nodes);

                for(const auto& qp : quadrature_points) {
                    Eigen::VectorXd dN_dxi = Utils::ShapeFunctions::getLineShapeFunctionDerivatives(order, qp.point(0));
                    Eigen::MatrixXd B = (1.0 / detJ) * dN_dxi.transpose();
                    ke_local += B.transpose() * local_sigma * B * qp.weight * detJ;
                }

                for(size_t i = 0; i < num_elem_nodes; ++i) {
                    for(size_t j = 0; j < num_elem_nodes; ++j) {
                        if (dofs[i] != -1 && dofs[j] != -1) {
                            k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i,j));
                        }
                    }
                }
            }
        }
        K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
        logger.info("Assembly for ", getName(), " complete.");
    }
}