#include <core/bcs/BoundaryCondition.hpp>
#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include <vector>
#include <limits>

namespace Core {

// --- DirichletBC Implementation (Direct Elimination Method) ---
DirichletBC::DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value, const std::string& tag)
    : BoundaryCondition(tag), value_(value) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}
DirichletBC::DirichletBC(int equation_index, Eigen::VectorXd value, const std::string& tag)
    : BoundaryCondition(tag), equation_index_(equation_index), value_(value) {
// This constructor directly accepts the equation index, no lookup needed.
}
void DirichletBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        Utils::Logger::instance().error("DirichletBC: Invalid equation index ", equation_index_);
        return;
    }

    // For all other rows j, update the RHS: F(j) = F(j) - K(j, i) * value
    // This must be done BEFORE zeroing out the column.
    for (int j = 0; j < K.rows(); ++j) {
        if (j == equation_index_) continue;
        F(j, 0) -= K.coeff(j, equation_index_) * value_(0);
    }

    // Zero out the row and column corresponding to the boundary condition.
    // This is a robust way to do it for a sparse matrix.
    for (int k=0; k<K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K,k); it; ++it)
        {
            if(it.row() == equation_index_ || it.col() == equation_index_) {
                it.valueRef() = 0;
            }
        }
    }

    // Set the diagonal element to 1 and the RHS to the prescribed value.
    K.coeffRef(equation_index_, equation_index_) = 1.0;
    F(equation_index_, 0) = value_(0);

    // It's good practice to prune the matrix of explicit zeros.
    K.prune(0.0);
}

// --- NeumannBC Implementation ---
NeumannBC::NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value, const std::string& tag)
        : BoundaryCondition(tag), flux_value_(flux_value) {
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

// --- CauchyBC Implementation ---
CauchyBC::CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf, const std::string& tag)
        : BoundaryCondition(tag), h_(h), T_inf_(T_inf) {
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


} // namespace Core
