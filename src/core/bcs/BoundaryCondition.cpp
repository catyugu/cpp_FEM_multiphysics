#include <core/bcs/BoundaryCondition.hpp>
#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include <vector>

namespace Core {

// --- DirichletBC Implementation ---
DirichletBC::DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value)
    : value_(value) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void DirichletBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        SimpleLogger::Logger::instance().error("DirichletBC: Invalid equation index ", equation_index_);
        return;
    }

    if (F.cols() != value_.size()) {
        SimpleLogger::Logger::instance().error("DirichletBC: Number of load cases in F does not match the size of the value vector.");
        return;
    }

    // Modify the RHS matrix for each load case
    for (int j = 0; j < K.outerSize(); ++j) {
        if (j == equation_index_) continue;
        double k_ji = K.coeff(j, equation_index_);
        if (k_ji != 0.0) {
            for (int col = 0; col < F.cols(); ++col) {
                F(j, col) -= k_ji * value_(col);
            }
        }
    }

    // Zero out the row and column in the stiffness matrix
    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            if (it.row() == equation_index_) it.valueRef() = 0.0;
        }
    }
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, equation_index_); it; ++it) {
        it.valueRef() = 0.0;
    }

    // Set the diagonal element to 1 and the RHS to the desired value
    K.coeffRef(equation_index_, equation_index_) = 1.0;
    F.row(equation_index_) = value_;
    K.prune(0.0);
}

// --- NeumannBC Implementation ---
NeumannBC::NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value)
    : flux_value_(flux_value) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void NeumannBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= F.rows()) {
        SimpleLogger::Logger::instance().error("NeumannBC: Invalid equation index ", equation_index_);
        return;
    }
    if (F.cols() != flux_value_.size()) {
        SimpleLogger::Logger::instance().error("NeumannBC: Number of load cases in F does not match the size of the flux_value vector.");
        return;
    }

    for (int i = 0; i < F.cols(); ++i) {
        F(equation_index_, i) += flux_value_(i);
    }
}

// --- CauchyBC Implementation ---
CauchyBC::CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf)
    : h_(h), T_inf_(T_inf) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void CauchyBC::apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        SimpleLogger::Logger::instance().error("CauchyBC: Invalid equation index ", equation_index_);
        return;
    }
    if (F.cols() != h_.size() || F.cols() != T_inf_.size()) {
        SimpleLogger::Logger::instance().error("CauchyBC: Number of load cases in F does not match the size of h or T_inf vectors.");
        return;
    }

    // The modification to the stiffness matrix K is often independent of the load case,
    // assuming 'h' is constant across them. If 'h' can vary per load case, this would need
    // a more complex solver structure. Here, we assume an average or primary 'h' value.
    K.coeffRef(equation_index_, equation_index_) += h_(0); // Using the first value as representative

    // Modify the RHS for each load case
    for (int i = 0; i < F.cols(); ++i) {
        F(equation_index_, i) += h_(i) * T_inf_(i);
    }
}


} // namespace Core