#include "solver/CoupledElectroThermalSolver.hpp"
#include "core/Problem.hpp"
#include "physics/Current3D.hpp"
#include "physics/Heat3D.hpp"
#include <solver/LinearSolver.hpp>
#include "utils/SimpleLogger.hpp"

namespace Solver {
    Eigen::VectorXd calculate_residual(const Eigen::SparseMatrix<double>& K, const Eigen::VectorXd& U, const Eigen::VectorXd& F) {
        return K * U - F;
    }


    void CoupledElectroThermalSolver::solveSteadyState(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem with Damped Newton-like Iterations ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();
        const auto& dof_manager = problem.getDofManager();
        const auto& nodes = problem.getMesh().getNodes();

        const double damping_factor = 0.8;

        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->getSolution().setConstant(293.15);
        }

        Eigen::VectorXd T_prev_iter = heat_field->getSolution();

        for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
            logger.info("--> Iteration ", iter + 1, " / ", problem.getMaxIterations());

            // --- Step 1: Solve EMag Field ---
            logger.info("    Solving EMag Field...");
            emag_field->assemble(heat_field);
            emag_field->applySources();

            Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
            Eigen::VectorXd F_emag_solve = emag_field->getRHS();

            for (const auto& node : nodes) {
                int temp_dof = dof_manager.getEquationIndex(node->getId(), "Temperature");
                if (temp_dof != -1) {
                    K_emag_solve.coeffRef(temp_dof, temp_dof) = 1.0;
                    F_emag_solve(temp_dof) = heat_field->getSolution()(temp_dof);
                }
            }

            for (const auto &bc: emag_field->getBCs()) {
                bc->apply(K_emag_solve, F_emag_solve);
            }
            LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

            // --- Step 2: Execute Coupling ---
            coupling_manager.executeCouplings();

            // --- Step 3: Solve Heat Field with a Corrected Newton Step ---
            logger.info("    Solving Heat Field (Newton-like step)...");
            heat_field->assemble();
            heat_field->applySources();

            Eigen::SparseMatrix<double> K_heat_tangent = heat_field->getStiffnessMatrix();
            Eigen::VectorXd F_heat_total = heat_field->getRHS();

            // **FIX START**: Correct handling of residual and boundary conditions
            // 1. Calculate the residual R = K(T_k) * U_k - F_total
            Eigen::VectorXd R = K_heat_tangent * heat_field->getSolution() - F_heat_total;

            // 2. Apply BCs to the tangent matrix (this modifies K_heat_tangent in place)
            //    We pass a dummy vector for the RHS because we'll handle the residual separately.
            Eigen::VectorXd dummy_rhs = Eigen::VectorXd::Zero(R.size());
            for (const auto &bc: heat_field->getBCs()) {
                bc->apply(K_heat_tangent, dummy_rhs);
            }

            // 3. Manually zero out the residual at the Dirichlet nodes. This is the correct approach.
            for (const auto &bc: heat_field->getBCs()) {
                if (auto* dirichlet = dynamic_cast<const Core::DirichletBC*>(bc.get())) {
                    int eq_index = dirichlet->getEquationIndex();
                    if(eq_index != -1) R(eq_index) = 0.0;
                }
            }

            // 4. Stabilize the matrix and residual for inactive DOFs
            for (const auto& node : nodes) {
                int volt_dof = dof_manager.getEquationIndex(node->getId(), "Voltage");
                if (volt_dof != -1) {
                    K_heat_tangent.coeffRef(volt_dof, volt_dof) = 1.0;
                    R(volt_dof) = 0.0;
                }
            }

            // 5. Solve for the update using the modified system K_t * delta_T = -R
            Eigen::VectorXd delta_T(heat_field->getSolution().size());
            LinearSolver::solve(K_heat_tangent, -R, delta_T);
            // **FIX END**

            heat_field->getSolution() += damping_factor * delta_T;

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
        logger.warn("--- Coupled steady-state solver did not converge after ", problem.getMaxIterations(), " iterations. ---");
    }

    void CoupledElectroThermalSolver::solveTransient(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Starting Coupled Transient Solve ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();
        const auto& dof_manager = problem.getDofManager();
        const auto& nodes = problem.getMesh().getNodes();

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

                // Stabilize for Temperature DOFs
                for (const auto& node : nodes) {
                    int temp_dof = dof_manager.getEquationIndex(node->getId(), "Temperature");
                    if (temp_dof != -1) {
                        K_emag_solve.coeffRef(temp_dof, temp_dof) = 1.0;
                        F_emag_solve(temp_dof) = heat_field->getSolution()(temp_dof);
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

                // Stabilize for Voltage DOFs
                for (const auto& node : nodes) {
                    int volt_dof = dof_manager.getEquationIndex(node->getId(), "Voltage");
                    if (volt_dof != -1) {
                        A_eff.coeffRef(volt_dof, volt_dof) = 1.0;
                        b_eff(volt_dof) = emag_field->getSolution()(volt_dof);
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