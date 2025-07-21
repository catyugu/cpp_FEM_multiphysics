#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers> // For BiCGSTAB
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"

namespace Solver {

    enum class SolverType {
        LU,
        BiCGSTAB
    };

    class LinearSolver {
    public:
        /**
         * @brief Solves the linear system Ax = b for x.
         * @param A The stiffness matrix.
         * @param b The right-hand side vector.
         * @param x The solution vector (output).
         * @param solver_type The type of solver to use (LU or BiCGSTAB).
         * @param max_iterations Maximum number of iterations for iterative solvers.
         * @param tolerance Convergence tolerance for iterative solvers.
         */
        static void solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x,
                          SolverType solver_type = SolverType::LU,
                          int max_iterations = 1000, double tolerance = 1e-9) {
            auto& logger = Utils::Logger::instance();
            logger.info("Starting linear solve...");

            if (solver_type == SolverType::LU) {
                Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
                solver.compute(A);
                if(solver.info() != Eigen::Success) {
                    throw Exception::SolverException("LU decomposition failed.");
                }

                x = solver.solve(b);
                if(solver.info() != Eigen::Success) {
                    throw Exception::SolverException("The solve step failed.");
                }
            } else if (solver_type == SolverType::BiCGSTAB) {
                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
                solver.setMaxIterations(max_iterations);
                solver.setTolerance(tolerance);
                solver.compute(A);
                if(solver.info() != Eigen::Success) {
                    throw Exception::SolverException("BiCGSTAB decomposition failed (preconditioner setup).");
                }

                x = solver.solve(b);
                if(solver.info() != Eigen::Success) {
                    throw Exception::SolverException("BiCGSTAB solve failed or did not converge.");
                }
                logger.info("BiCGSTAB converged in ", solver.iterations(), " iterations with error ", solver.error());
            } else {
                throw Exception::ConfigurationException("Unknown solver type specified.");
            }

            logger.info("Linear solve successful.");
        }
    };

} // namespace Solver

#endif // LINEARSOLVER_HPP