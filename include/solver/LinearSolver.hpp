// include/solver/LinearSolver.hpp (修改后的版本)

#ifndef LINEARSOLVER_HPP
#define LINEARSOLVER_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers> // For BiCGSTAB
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"
#include <chrono> // 用于计时

namespace Solver {

    enum class SolverType {
        LU,
        BiCGSTAB
    };

    class LinearSolver {
    public:
        static void solve(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x,
                          SolverType solver_type = SolverType::LU,
                          int max_iterations = 1000, double tolerance = 1e-9) {
            auto& logger = Utils::Logger::instance();
            logger.info("Starting linear solve...");
            auto start_time = std::chrono::high_resolution_clock::now();

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
                // --- 核心修改：使用 IncompleteLUT 预条件器 ---
                Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
                solver.setMaxIterations(max_iterations);
                solver.setTolerance(tolerance);

                logger.info("    BiCGSTAB: Computing preconditioner...");
                solver.compute(A); // compute() 会计算预条件器
                if(solver.info() != Eigen::Success) {
                    throw Exception::SolverException("BiCGSTAB compute/preconditioner setup failed.");
                }
                logger.info("    BiCGSTAB: Solving...");
                x = solver.solve(b);

                if(solver.info() != Eigen::Success) {
                    throw Exception::SolverException("BiCGSTAB solve failed or did not converge.");
                }
                logger.info("BiCGSTAB converged in ", solver.iterations(), " iterations with error ", solver.error());
            } else {
                throw Exception::ConfigurationException("Unknown solver type specified.");
            }

            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            logger.info("Linear solve successful. Time taken: ", elapsed.count(), "s");
        }
    };

} // namespace Solver

#endif // LINEARSOLVER_HPP