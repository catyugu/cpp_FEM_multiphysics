#include <core/bcs/BoundaryCondition.hpp>
#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include <vector>
#include <limits>
#include "core/mesh/Mesh.hpp"
#include "core/mesh/Element.hpp"
#include <set>
#include <algorithm>
#include <cmath>

namespace Core {

// --- DirichletBC Implementation (Restored to full functionality) ---
DirichletBC::DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value, const std::string& tag)
    : BoundaryCondition(tag), value_(std::move(value)) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

DirichletBC::DirichletBC(int equation_index, Eigen::VectorXd value, const std::string& tag)
    : BoundaryCondition(tag), equation_index_(equation_index), value_(std::move(value)) {}

void DirichletBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        Utils::Logger::instance().error("DirichletBC: Invalid equation index ", equation_index_);
        return;
    }

    // 1. 修改力向量 F
    for (int j = 0; j < K.rows(); ++j) {
        if (j != equation_index_) {
            F(j, 0) -= K.coeff(j, equation_index_) * value_(0);
        }
    }

    // 2. 将K矩阵的对应行和列清零
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            if (it.row() == equation_index_ || it.col() == equation_index_) {
                it.valueRef() = 0.0;
            }
        }
    }

    // 3. 设置对角线元素为1，并设置F向量的对应值
    K.coeffRef(equation_index_, equation_index_) = 1.0;
    F(equation_index_, 0) = value_(0);

    // 移除零元素以保持矩阵的稀疏性
    K.prune(0.0);
}


// --- NeumannBC and CauchyBC Implementations (Unchanged) ... ---
NeumannBC::NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value, const std::string& tag)
        : BoundaryCondition(tag), flux_value_(std::move(flux_value)) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void NeumannBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= F.rows()) {
        Utils::Logger::instance().error("NeumannBC: Invalid equation index ", equation_index_);
        return;
    }
    F(equation_index_, 0) += flux_value_(0);
}

CauchyBC::CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf, const std::string& tag)
        : BoundaryCondition(tag), h_(std::move(h)), T_inf_(std::move(T_inf)) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void CauchyBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        Utils::Logger::instance().error("CauchyBC: Invalid equation index ", equation_index_);
        return;
    }
    K.coeffRef(equation_index_, equation_index_) += h_(0);
    F(equation_index_, 0) += h_(0) * T_inf_(0);
}


// --- Factory Method Implementation (Unchanged) ---
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

                    if (region_predicate(node1->getCoords()) && region_predicate(node2->getCoords())) {
                         bool edge_is_truly_on_boundary = false;
                         std::vector<double> midpoint_coords = {
                            (node1->getCoords()[0] + node2->getCoords()[0]) / 2.0,
                            (node1->getCoords()[1] + node2->getCoords()[1]) / 2.0,
                            (node1->getCoords().size() > 2 ? (node1->getCoords()[2] + node2->getCoords()[2]) / 2.0 : 0.0)
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