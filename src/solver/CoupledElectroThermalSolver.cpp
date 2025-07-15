#include "solver/CoupledElectroThermalSolver.hpp"
#include "core/Problem.hpp"
#include "physics/Current1D.hpp"
#include "physics/Heat1D.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include <solver/LinearSolver.hpp>
#include "utils/SimpleLogger.hpp"

namespace Solver {
    void CoupledElectroThermalSolver::solveSteadyState(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();
        const auto &dof_manager = problem.getDofManager();

        // Initialize temperature solution to a reasonable default if it's zero
        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->getSolution().setConstant(293.15); // Room temperature
        }
        Eigen::MatrixXd T_prev_iter = heat_field->getSolution();

        for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
            logger.info("--> Iteration ", iter + 1, " / ", problem.getMaxIterations());

            // --- Step 1: Solve EMag Field ---
            logger.info("    Solving EMag Field...");
            emag_field->assemble();
            emag_field->applySources();

            // Create temporary copies to build a well-posed system for this step
            Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
            Eigen::MatrixXd F_emag_solve = emag_field->getRHS();

            // Constrain Temperature DOFs for this EMag solve step
            // This ensures the overall system matrix is not singular
            for (const auto &node: problem.getMesh().getNodes()) {
                int dof_idx = dof_manager.getEquationIndex(node->getId(), "Temperature");
                if (dof_idx != -1) {
                    K_emag_solve.coeffRef(dof_idx, dof_idx) = 1.0;
                    F_emag_solve(dof_idx, 0) = heat_field->getSolution()(dof_idx, 0);
                }
            }

            // Apply EMag BCs to the temporary system
            for (const auto &bc: emag_field->getBCs()) {
                bc->apply(K_emag_solve, F_emag_solve);
            }

            // Solve the temporary system, updating the field's official solution vector
            LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

            // --- Step 2: Execute Coupling ---
            coupling_manager.executeCouplings();

            // --- Step 3: Solve Heat Field ---
            logger.info("    Solving Heat Field...");
            heat_field->assemble();
            heat_field->applySources(); // Includes Joule heating

            // Create temporary copies for the heat solve
            Eigen::SparseMatrix<double> K_heat_solve = heat_field->getStiffnessMatrix();
            Eigen::MatrixXd F_heat_solve = heat_field->getRHS();

            // Constrain Voltage DOFs for this Heat solve step
            for (const auto &node: problem.getMesh().getNodes()) {
                int dof_idx = dof_manager.getEquationIndex(node->getId(), "Voltage");
                if (dof_idx != -1) {
                    K_heat_solve.coeffRef(dof_idx, dof_idx) = 1.0;
                    F_heat_solve(dof_idx, 0) = emag_field->getSolution()(dof_idx, 0);
                }
            }

            // Apply Heat BCs to the temporary system
            for (const auto &bc: heat_field->getBCs()) {
                bc->apply(K_heat_solve, F_heat_solve);
            }

            LinearSolver::solve(K_heat_solve, F_heat_solve, heat_field->getSolution());

            // --- Step 4: Check for Convergence ---
            double norm_diff = (heat_field->getSolution() - T_prev_iter).norm();
            double norm_sol = heat_field->getSolution().norm();
            double relative_error = (norm_sol > 1e-12) ? (norm_diff / norm_sol) : norm_diff;

            logger.info("    Convergence Check: Relative Error = ", relative_error);
            if (relative_error < problem.getConvergenceTolerance()) {
                logger.info("--- Coupled steady-state solver converged after ", iter + 1, " iterations. ---");
                return;
            }

            T_prev_iter = heat_field->getSolution();
        }

        logger.warn("--- Coupled steady-state solver did not converge after ", problem.getMaxIterations(),
                    " iterations. ---");
    }


    void CoupledElectroThermalSolver::solveTransient(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Starting Coupled Transient Solve (Mixed-Order Robust) ---");
        logger.info("Time Step: ", problem.getTimeStep(), "s, Total Time: ", problem.getTotalTime(), "s");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();
        const auto &dof_manager = problem.getDofManager();
        const auto &all_vars = dof_manager.getVariableNames();

        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->setInitialConditions(293.15);
        }

        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        double dt = problem.getTimeStep();

        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * dt, "s");

            heat_field->updatePreviousSolution();
            Eigen::MatrixXd T_prev_iter_inner_loop = heat_field->getSolution();

            for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
                logger.info("  --> Inner Iteration ", iter + 1, " / ", problem.getMaxIterations());

                // --- Step 1: Solve EMag Field ---
                emag_field->assemble();
                emag_field->applySources();

                Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
                Eigen::MatrixXd F_emag_solve = emag_field->getRHS();

                // Stabilize all DOFs that are NOT for "Voltage"
                for (const auto &var_name: all_vars) {
                    if (var_name != emag_field->getVariableName()) {
                        for (int dof_idx = 0; dof_idx < dof_manager.getNumEquations(); ++dof_idx) {
                            // A simple stabilization for now
                            K_emag_solve.coeffRef(dof_idx, dof_idx) += 1.0;
                        }
                    }
                }
                emag_field->applyBCs(); // Apply EMag BCs after stabilization

                LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

                // --- Step 2: Execute Coupling ---
                coupling_manager.executeCouplings();

                // --- Step 3: Solve Heat Field ---
                heat_field->assemble();
                heat_field->applySources();

                Eigen::SparseMatrix<double> A_eff =
                        (heat_field->getMassMatrix() / dt) + heat_field->getStiffnessMatrix();
                Eigen::MatrixXd b_eff = heat_field->getRHS() + (heat_field->getMassMatrix() / dt) * heat_field->
                                        getPreviousSolution();

                // Stabilize all DOFs that are NOT for "Temperature"
                for (const auto &var_name: all_vars) {
                    if (var_name != heat_field->getVariableName()) {
                        for (int dof_idx = 0; dof_idx < dof_manager.getNumEquations(); ++dof_idx) {
                            A_eff.coeffRef(dof_idx, dof_idx) += 1.0;
                        }
                    }
                }
                heat_field->applyBCs(); // Apply Heat BCs after stabilization

                LinearSolver::solve(A_eff, b_eff, heat_field->getSolution());

                // --- Step 4: Check for Inner Loop Convergence ---
                double norm_diff = (heat_field->getSolution() - T_prev_iter_inner_loop).norm();
                double norm_sol = heat_field->getSolution().norm();
                double relative_error = (norm_sol > 1e-12) ? (norm_diff / norm_sol) : norm_diff;

                logger.info("    Inner Loop Convergence Check: Relative Error = ", relative_error);
                if (relative_error < problem.getConvergenceTolerance()) {
                    logger.info("  Inner loop converged for time step ", i + 1, " after ", iter + 1, " iterations.");
                    break;
                }
                T_prev_iter_inner_loop = heat_field->getSolution();
            }
        }
        logger.info("\n--- Coupled Transient Solve Finished ---");
    }
} // Solver
