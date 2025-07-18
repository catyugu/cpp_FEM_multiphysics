#include <core/bcs/BoundaryCondition.hpp>
#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include <vector>
#include <limits>
#include "core/mesh/Mesh.hpp" // Added include
#include "core/mesh/Element.hpp" // Added include
#include <set> // Added include
#include <algorithm> // Added include
#include <cmath> // For std::abs

namespace Core {

// --- DirichletBC Implementation (Direct Elimination Method) ---
DirichletBC::DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value, const std::string& tag)
    : BoundaryCondition(tag), value_(std::move(value)) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}
DirichletBC::DirichletBC(int equation_index, Eigen::VectorXd value, const std::string& tag)
    : BoundaryCondition(tag), equation_index_(equation_index), value_(std::move(value)) {
// This constructor directly accepts the equation index, no lookup needed.
}
void DirichletBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        Utils::Logger::instance().error("DirichletBC: Invalid equation index ", equation_index_);
        return;
    }

    for (int j = 0; j < K.rows(); ++j) {
        if (j == equation_index_) continue;
        F(j, 0) -= K.coeff(j, equation_index_) * value_(0);
    }

    for (int k=0; k<K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K,k); it; ++it)
        {
            if(it.row() == equation_index_ || it.col() == equation_index_) {
                it.valueRef() = 0;
            }
        }
    }

    K.coeffRef(equation_index_, equation_index_) = 1.0;
    F(equation_index_, 0) = value_(0);
    K.prune(0.0);
}

// --- NeumannBC and CauchyBC Implementations ... ---
NeumannBC::NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value, const std::string& tag)
        : BoundaryCondition(tag), flux_value_(std::move(flux_value)) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void NeumannBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= F.rows()) {
        Utils::Logger::instance().error("NeumannBC: Invalid equation index ", equation_index_);
        return;
    }
    if (F.cols() != flux_value_.size()) {
        Utils::Logger::instance().error("NeumannBC: Number of load cases in F does not match the size of the flux_value vector.");
        return;
    }

    for (int i = 0; i < F.cols(); ++i) {
        F(equation_index_, i) += flux_value_(i);
    }
}

CauchyBC::CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf, const std::string& tag)
        : BoundaryCondition(tag), h_(std::move(h)), T_inf_(std::move(T_inf)) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void CauchyBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        Utils::Logger::instance().error("CauchyBC: Invalid equation index ", equation_index_);
        return;
    }
    if (F.cols() != h_.size() || F.cols() != T_inf_.size()) {
        Utils::Logger::instance().error("CauchyBC: Number of load cases in F does not match the size of h or T_inf vectors.");
        return;
    }

    K.coeffRef(equation_index_, equation_index_) += h_(0);
    for (int i = 0; i < F.cols(); ++i) {
        F(equation_index_, i) += h_(i) * T_inf_(i);
    }
}


// --- New Factory Method Implementation ---
std::vector<std::unique_ptr<BoundaryCondition>> DirichletBC::create(
    const DOFManager& dof_manager,
    const Mesh& mesh,
    const std::string& var_name,
    int field_order,
    const std::function<bool(const std::vector<double>& coords)>& region_predicate,
    const std::function<double(const std::vector<double>& coords)>& value_function,
    const std::string& tag
) {
    auto& logger = Utils::Logger::instance();
    logger.info("Creating DirichletBCs from predicate for variable '", var_name, "' with order ", field_order);

    std::vector<std::unique_ptr<BoundaryCondition>> bcs;
    std::set<int> constrained_dofs;

    for (const auto& node : mesh.getNodes()) {
        if (region_predicate(node->getCoords())) {
            int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
            if (dof_idx != -1 && constrained_dofs.find(dof_idx) == constrained_dofs.end()) {
                double bc_value = value_function(node->getCoords());
                bcs.push_back(std::make_unique<DirichletBC>(dof_idx, Eigen::Vector<double, 1>(bc_value), tag));
                constrained_dofs.insert(dof_idx);
            }
        }
    }

    if (field_order > 1) {
        for (const auto& elem : mesh.getElements()) {
            const auto& elem_nodes = elem->getNodes();
            for (size_t i = 0; i < elem_nodes.size(); ++i) {
                for (size_t j = i + 1; j < elem_nodes.size(); ++j) {
                    auto* node1 = elem_nodes[i];
                    auto* node2 = elem_nodes[j];

                    // *** THIS IS THE FIX ***
                    // An edge is on the boundary if both its vertices satisfy the predicate AND
                    // the edge lies flat on one of the primary boundary planes (for axis-aligned geometry).
                    if (region_predicate(node1->getCoords()) && region_predicate(node2->getCoords())) {
                         bool edge_is_truly_on_boundary = false;
                         constexpr double coord_eps = 1e-9;
                         // Check if nodes share a common x, y, or z coordinate that is on the boundary
                         if (std::abs(node1->getCoords()[0] - node2->getCoords()[0]) < coord_eps && region_predicate(node1->getCoords())) {
                            edge_is_truly_on_boundary = true;
                         } else if (std::abs(node1->getCoords()[1] - node2->getCoords()[1]) < coord_eps && region_predicate(node1->getCoords())) {
                            edge_is_truly_on_boundary = true;
                         } else if (elem->getDimension() == 3 && std::abs(node1->getCoords()[2] - node2->getCoords()[2]) < coord_eps && region_predicate(node1->getCoords())) {
                            edge_is_truly_on_boundary = true;
                         }

                        // A more general (but slower) check: is the midpoint also on the boundary?
                        std::vector<double> midpoint_coords = {
                            (node1->getCoords()[0] + node2->getCoords()[0]) / 2.0,
                            (node1->getCoords()[1] + node2->getCoords()[1]) / 2.0,
                            (node1->getCoords()[2] + node2->getCoords()[2]) / 2.0
                        };
                        if (region_predicate(midpoint_coords)) {
                             edge_is_truly_on_boundary = true;
                        }


                        if (edge_is_truly_on_boundary) {
                            std::vector<int> edge_node_ids = {node1->getId(), node2->getId()};
                            int edge_dof_idx = dof_manager.getEdgeEquationIndex(edge_node_ids, var_name);

                            if (edge_dof_idx != -1 && constrained_dofs.find(edge_dof_idx) == constrained_dofs.end()) {
                                double bc_value = value_function(midpoint_coords);
                                bcs.push_back(std::make_unique<DirichletBC>(edge_dof_idx, Eigen::Vector<double, 1>(bc_value), tag));
                                constrained_dofs.insert(edge_dof_idx);
                            }
                        }
                    }
                }
            }
        }
    }

    logger.info("Generated ", bcs.size(), " DirichletBC objects from predicate.");
    return bcs;
}


} // namespace Core