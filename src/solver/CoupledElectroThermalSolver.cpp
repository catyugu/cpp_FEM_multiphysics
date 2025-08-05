#include "solver/CoupledElectroThermalSolver.hpp"
#include "core/Problem.hpp"
#include <solver/LinearSolver.hpp>
#include "utils/SimpleLogger.hpp"

namespace Solver {

    // ... (solveSteadyState logic remains the same as the last working version)
    void CoupledElectroThermalSolver::solveSteadyState(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem with Damped Newton-like Iterations ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();
        const auto& mesh = problem.getMesh();
        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->getSolution().setConstant(293.15);
        }
        Eigen::VectorXd T_prev_iter = heat_field->getSolution();

        for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
            constexpr double damping_factor = 0.8;
            logger.info("--> Iteration ", iter + 1, " / ", problem.getMaxIterations());
            logger.info("    Solving EMag Field...");
            emag_field->assemble(heat_field);
            emag_field->applySources();
            Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
            Eigen::VectorXd F_emag_solve = emag_field->getRHS();
            for (const auto& elem : mesh.getElements()) {
                const auto heat_element_dofs = heat_field->getElementDofs(elem);
                for (int dof_idx : heat_element_dofs) {
                    if (dof_idx != -1) {
                        K_emag_solve.coeffRef(dof_idx, dof_idx) = 1.0;
                        F_emag_solve(dof_idx) = heat_field->getSolution()(dof_idx);
                    }
                }
            }
            for (const auto &bc: emag_field->getBCs()) {
                bc->apply(K_emag_solve, F_emag_solve);
            }
            LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

            coupling_manager.executeCouplings();

            logger.info("    Solving Heat Field (Newton-like step)...");
            heat_field->assemble();
            heat_field->applySources();

            Eigen::SparseMatrix<double> K_heat_tangent = heat_field->getStiffnessMatrix();
            Eigen::VectorXd F_heat_total = heat_field->getRHS() + heat_field->getCouplingRHS();
            Eigen::VectorXd R = K_heat_tangent * heat_field->getSolution() - F_heat_total;

            Eigen::VectorXd dummy_rhs = Eigen::VectorXd::Zero(R.size());
            for (const auto &bc: heat_field->getBCs()) {
                bc->apply(K_heat_tangent, dummy_rhs);
            }
            for (const auto &bc: heat_field->getBCs()) {
                if (auto* dirichlet = dynamic_cast<const Core::DirichletBC*>(bc.get())) {
                    int eq_index = dirichlet->getEquationIndex();
                    if(eq_index != -1) R(eq_index) = 0.0;
                }
            }

            for (const auto& elem : mesh.getElements()) {
                 const auto emag_element_dofs = emag_field->getElementDofs(elem);
                 for (int dof_idx : emag_element_dofs) {
                    if (dof_idx != -1) {
                        K_heat_tangent.coeffRef(dof_idx, dof_idx) = 1.0;
                        R(dof_idx) = 0.0;
                    }
                 }
            }

            Eigen::VectorXd delta_T(heat_field->getSolution().size());
            LinearSolver::solve(K_heat_tangent, -R, delta_T);
            heat_field->getSolution() += damping_factor * delta_T;

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
        const auto& mesh = problem.getMesh();
        if (heat_field->getSolution().isZero(1e-9)) {
            heat_field->setInitialConditions(293.15);
        }
        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        double dt = problem.getTimeStep();
        emag_field->updatePreviousSolution();
        heat_field->updatePreviousSolution();

        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * dt, "s");
            Eigen::VectorXd T_prev_iter_inner = heat_field->getSolution();

            for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
                logger.info("  --> Inner Iteration ", iter + 1, " / ", problem.getMaxIterations());

                emag_field->assemble(heat_field);
                emag_field->applySources();
                Eigen::SparseMatrix<double> K_emag_solve = emag_field->getStiffnessMatrix();
                Eigen::VectorXd F_emag_solve = emag_field->getRHS();
                for (const auto& elem : mesh.getElements()) {
                    for (int dof_idx : heat_field->getElementDofs(elem)) {
                        if (dof_idx != -1) {
                            K_emag_solve.coeffRef(dof_idx, dof_idx) = 1.0;
                            F_emag_solve(dof_idx) = heat_field->getSolution()(dof_idx);
                        }
                    }
                }
                for (const auto& bc : emag_field->getBCs()) {
                    bc->apply(K_emag_solve, F_emag_solve);
                }
                LinearSolver::solve(K_emag_solve, F_emag_solve, emag_field->getSolution());

                coupling_manager.executeCouplings();

                heat_field->assemble();
                heat_field->applySources();

                Eigen::SparseMatrix<double> A_eff = (heat_field->getMassMatrix() / dt) + heat_field->getStiffnessMatrix();
                Eigen::VectorXd b_eff = heat_field->getRHS() + heat_field->getCouplingRHS() + (heat_field->getMassMatrix() / dt) * heat_field->getPreviousSolution();

                for (const auto& elem : mesh.getElements()) {
                    for (int dof_idx : emag_field->getElementDofs(elem)) {
                        if (dof_idx != -1) {
                            A_eff.coeffRef(dof_idx, dof_idx) = 1.0;
                            b_eff(dof_idx) = emag_field->getSolution()(dof_idx);
                        }
                    }
                }
                for (const auto& bc : heat_field->getBCs()) {
                    bc->apply(A_eff, b_eff);
                }
                LinearSolver::solve(A_eff, b_eff, heat_field->getSolution());

                double norm_diff = (heat_field->getSolution() - T_prev_iter_inner).norm();
                double norm_sol = heat_field->getSolution().norm();
                double relative_error = (norm_sol > 1e-12) ? (norm_diff / norm_sol) : norm_diff;
                logger.info("    Inner Loop Convergence Check: Relative Error = ", relative_error);
                if (relative_error < problem.getConvergenceTolerance()) {
                    logger.info("  Inner loop converged for time step ", i + 1, " after ", iter + 1, " iterations.");
                    break;
                }
                 if (iter == problem.getMaxIterations() - 1) {
                    logger.warn("  Inner loop did not converge for time step ", i + 1);
                }
                T_prev_iter_inner = heat_field->getSolution();
            }

            heat_field->updatePreviousSolution();
            emag_field->updatePreviousSolution();
        }
        logger.info("\n--- Coupled Transient Solve Finished ---");
    }

} // namespace Solver