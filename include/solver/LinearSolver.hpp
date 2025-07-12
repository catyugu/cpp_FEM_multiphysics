#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"

namespace Solver {

    class LinearSolver {
    public:
        // Solves the system Ax = b and returns x
        static void solve(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& b, Eigen::MatrixXd& x) {
            auto& logger = SimpleLogger::Logger::instance();
            logger.info("Starting linear solve...");

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            if(solver.info() != Eigen::Success) {
                throw Exception::SolverException("LU decomposition failed.");
            }

            x = solver.solve(b);
            if(solver.info() != Eigen::Success) {
                throw Exception::SolverException("The solve step failed.");
            }

            logger.info("Linear solve successful.");
        }
    };

} // namespace Solver


#endif // LINEARSOLVER_HPP