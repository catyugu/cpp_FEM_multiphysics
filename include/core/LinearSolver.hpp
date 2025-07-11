#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "utils/SimpleLogger.hpp"

namespace Core {

    class LinearSolver {
    public:
        // Solves the system Ax = b and returns x
        static bool solve(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& b, Eigen::MatrixXd& x) {
            auto& logger = SimpleLogger::Logger::instance();
            logger.info("Starting linear solve...");

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            if(solver.info() != Eigen::Success) {
                logger.error("Linear solver: LU decomposition failed.");
                return false;
            }

            x = solver.solve(b);
            if(solver.info() != Eigen::Success) {
                logger.error("Linear solver: The solve step failed.");
                return false;
            }

            logger.info("Linear solve successful.");
            return true;
        }
    };

} // namespace Core


#endif // LINEARSOLVER_HPP
