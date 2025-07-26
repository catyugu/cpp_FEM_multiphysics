#include "solver/CoupledElectroThermalSolver.hpp"
#include "core/Problem.hpp"
#include "physics/Current3D.hpp"
#include "physics/Heat3D.hpp"
#include <solver/LinearSolver.hpp>
#include "utils/SimpleLogger.hpp"

namespace Solver {
    void CoupledElectroThermalSolver::solveSteadyState(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();

        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->getSolution().setConstant(293.15);
        }
        Eigen::MatrixXd T_prev_iter = heat_field->getSolution();

        for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
            logger.info("--> Iteration ", iter + 1, " / ", problem.getMaxIterations());

            // --- Step 1: Solve EMag Field ---
            logger.info("    Solving EMag Field...");
            // FIX: Pass the heat_field to assemble() to provide temperature context
            emag_field->assemble(heat_field);
            emag_field->applySources();

            Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
            Eigen::VectorXd F_emag_solve = emag_field->getRHS();

            for (const auto& elem : problem.getMesh().getElements()) {
                elem->setOrder(heat_field->getElementOrder());
                // FIX: Corrected method name to get_element_dofs
                const auto heat_element_dofs = heat_field->getElementDofs(elem);

                for (int dof_idx_local : heat_element_dofs) {
                    if (dof_idx_local != -1) {
                        K_emag_solve.coeffRef(dof_idx_local, dof_idx_local) = 1.0;
                        F_emag_solve(dof_idx_local, 0) = heat_field->getSolution()(dof_idx_local, 0);
                    }
                }
            }

            for (const auto &bc: emag_field->getBCs()) {
                bc->apply(K_emag_solve, F_emag_solve);
            }

            LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

            // --- Step 2: Execute Coupling ---
            coupling_manager.executeCouplings();

            // --- Step 3: Solve Heat Field ---
            logger.info("    Solving Heat Field...");
            heat_field->assemble();
            heat_field->applySources();

            Eigen::SparseMatrix<double> K_heat_solve = heat_field->getStiffnessMatrix();
            Eigen::VectorXd F_heat_solve = heat_field->getRHS();

            for (const auto& elem : problem.getMesh().getElements()) {
                elem->setOrder(emag_field->getElementOrder());
                 // FIX: Corrected method name to get_element_dofs
                const auto emag_element_dofs = emag_field->getElementDofs(elem);

                for (int dof_idx_local : emag_element_dofs) {
                    if (dof_idx_local != -1) {
                        K_heat_solve.coeffRef(dof_idx_local, dof_idx_local) = 1.0;
                        F_heat_solve(dof_idx_local, 0) = emag_field->getSolution()(dof_idx_local, 0);
                    }
                }
            }

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
                // FIX: Removed the line that overwrites the voltage solution
                return;
            }
            T_prev_iter = heat_field->getSolution();
        }
        logger.warn("--- Coupled steady-state solver did not converge after ", problem.getMaxIterations(), " iterations. ---");
    }

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
                emag_field->assemble(heat_field);
                emag_field->applySources();
                Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
                Eigen::VectorXd F_emag_solve = emag_field->getRHS();

                for (const auto& elem : problem.getMesh().getElements()) {
                    elem->setOrder(heat_field->getElementOrder());
                    const auto heat_element_dofs = heat_field->getElementDofs(elem);
                    for (int dof_idx_local : heat_element_dofs) {
                        if (dof_idx_local != -1) {
                            K_emag_solve.coeffRef(dof_idx_local, dof_idx_local) = 1.0;
                            F_emag_solve(dof_idx_local, 0) = heat_field->getSolution()(dof_idx_local, 0);
                        }
                    }
                }

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
                    const auto emag_element_dofs = emag_field->getElementDofs(elem);
                    for (int dof_idx_local : emag_element_dofs) {
                        if (dof_idx_local != -1) {
                            A_eff.coeffRef(dof_idx_local, dof_idx_local) = 1.0;
                            b_eff(dof_idx_local, 0) = emag_field->getSolution()(dof_idx_local, 0);
                        }
                    }
                }

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