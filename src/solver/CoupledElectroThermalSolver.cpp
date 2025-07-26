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
        // const auto &all_vars = dof_manager.getVariableNames(); // No longer needed for this stabilization approach

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
            Eigen::VectorXd F_emag_solve = emag_field->getRHS();

            // FIX: Stabilize Temperature DOFs using their current solution values
            // This ensures the overall system matrix is not singular and respects current state.
            for (const auto& elem : problem.getMesh().getElements()) {
                // Get DOFs for the Temperature field on this element, respecting its order
                elem->setOrder(heat_field->getElementOrder());
                const auto heat_element_dofs = heat_field->get_element_dofs(elem);

                for (int dof_idx_local : heat_element_dofs) {
                    if (dof_idx_local != -1) {
                        K_emag_solve.coeffRef(dof_idx_local, dof_idx_local) = 1.0; // Set diagonal to 1.0
                        // Set RHS to current temperature solution value
                        F_emag_solve(dof_idx_local, 0) = heat_field->getSolution()(dof_idx_local, 0);
                    }
                }
            }

            // Apply EMag BCs to the temporary system (after stabilization)
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
            Eigen::VectorXd F_heat_solve = heat_field->getRHS();

            // FIX: Stabilize Voltage DOFs using their current solution values
            for (const auto& elem : problem.getMesh().getElements()) {
                // Get DOFs for the Voltage field on this element, respecting its order
                elem->setOrder(emag_field->getElementOrder());
                const auto emag_element_dofs = emag_field->get_element_dofs(elem);

                for (int dof_idx_local : emag_element_dofs) {
                    if (dof_idx_local != -1) {
                        K_heat_solve.coeffRef(dof_idx_local, dof_idx_local) = 1.0; // Set diagonal to 1.0
                        // Set RHS to current voltage solution value
                        F_heat_solve(dof_idx_local, 0) = emag_field->getSolution()(dof_idx_local, 0);
                    }
                }
            }

            // Apply Heat BCs to the temporary system (after stabilization)
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
                // **FIX**: Ensure the EMag field is updated with the final, consistent solution vector.
                problem.getField("Voltage")->getSolution() = problem.getField("Temperature")->getSolution();
                return;
            }
            T_prev_iter = heat_field->getSolution();
        }
        logger.warn("--- Coupled steady-state solver did not converge after ", problem.getMaxIterations(), " iterations. ---");
    }

    // The solveTransient method also needs the same stabilization fix for consistency
      void CoupledElectroThermalSolver::solveTransient(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Starting Coupled Transient Solve ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();

        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->setInitialConditions(293.15);
        }

        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        double dt = problem.getTimeStep();

        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * dt, "s");


            heat_field->updatePreviousSolution();
            emag_field->updatePreviousSolution();

            Eigen::MatrixXd T_prev_iter_inner_loop = heat_field->getSolution();

            for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
                logger.info("  --> Inner Iteration ", iter + 1, " / ", problem.getMaxIterations());

                // --- Step 1: Solve EMag Field ---
                emag_field->assemble();
                emag_field->applySources();
                Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
                Eigen::VectorXd F_emag_solve = emag_field->getRHS();

                for (const auto& elem : problem.getMesh().getElements()) {
                    elem->setOrder(heat_field->getElementOrder());
                    const auto heat_element_dofs = heat_field->get_element_dofs(elem);
                    for (int dof_idx_local : heat_element_dofs) {
                        if (dof_idx_local != -1) {
                            K_emag_solve.coeffRef(dof_idx_local, dof_idx_local) = 1.0;
                            F_emag_solve(dof_idx_local, 0) = heat_field->getSolution()(dof_idx_local, 0);
                        }
                    }
                }

                // **CORRECTION**: Apply BCs to the local matrices
                for (const auto& bc : emag_field->getBCs()) {
                    bc->apply(K_emag_solve, F_emag_solve);
                }

                LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

                // --- Step 2: Execute Coupling ---
                coupling_manager.executeCouplings();

                // --- Step 3: Solve Heat Field ---
                heat_field->assemble();
                heat_field->applySources();
                Eigen::SparseMatrix<double> A_eff = (heat_field->getMassMatrix() / dt) + heat_field->getStiffnessMatrix();
                Eigen::VectorXd b_eff = heat_field->getRHS() + (heat_field->getMassMatrix() / dt) * heat_field->getPreviousSolution();

                for (const auto& elem : problem.getMesh().getElements()) {
                    elem->setOrder(emag_field->getElementOrder());
                    const auto emag_element_dofs = emag_field->get_element_dofs(elem);
                    for (int dof_idx_local : emag_element_dofs) {
                        if (dof_idx_local != -1) {
                            A_eff.coeffRef(dof_idx_local, dof_idx_local) = 1.0;
                            b_eff(dof_idx_local, 0) = emag_field->getSolution()(dof_idx_local, 0);
                        }
                    }
                }

                // **CORRECTION**: Apply BCs to the local matrices
                for (const auto& bc : heat_field->getBCs()) {
                    bc->apply(A_eff, b_eff);
                }

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