#include "core/BoundaryCondition.hpp"
#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include <vector>

namespace Core {

// --- DirichletBC Implementation ---
DirichletBC::DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double value)
    : value_(value) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void DirichletBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        SimpleLogger::Logger::instance().error("DirichletBC: Invalid equation index ", equation_index_);
        return;
    }

    for (int j = 0; j < K.outerSize(); ++j) {
        if (j == equation_index_) continue;
        double k_ji = K.coeff(j, equation_index_);
        if (k_ji != 0.0) {
            F(j) -= k_ji * value_;
        }
    }

    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            if (it.row() == equation_index_) it.valueRef() = 0.0;
        }
    }
    for (Eigen::SparseMatrix<double>::InnerIterator it(K, equation_index_); it; ++it) {
        it.valueRef() = 0.0;
    }

    K.coeffRef(equation_index_, equation_index_) = 1.0;
    F(equation_index_) = value_;
    K.prune(0.0);
}

// --- NeumannBC Implementation ---
NeumannBC::NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double flux_value)
    : flux_value_(flux_value) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void NeumannBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= F.size()) {
        SimpleLogger::Logger::instance().error("NeumannBC: Invalid equation index ", equation_index_);
        return;
    }
    F(equation_index_) += flux_value_;
}

// --- CauchyBC Implementation ---
CauchyBC::CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double h, double T_inf)
    : h_(h), T_inf_(T_inf) {
    equation_index_ = dof_manager.getEquationIndex(node_id, var_name);
}

void CauchyBC::apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const {
    if (equation_index_ < 0 || equation_index_ >= K.rows()) {
        SimpleLogger::Logger::instance().error("CauchyBC: Invalid equation index ", equation_index_);
        return;
    }
    K.coeffRef(equation_index_, equation_index_) += h_;
    F(equation_index_) += h_ * T_inf_;
}


} // namespace Core
