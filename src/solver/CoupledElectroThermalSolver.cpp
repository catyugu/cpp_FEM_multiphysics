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
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();
        const auto &dof_manager = problem.getDofManager();
        const auto &mesh = problem.getMesh();

        Eigen::MatrixXd T_prev = heat_field->getSolution();
        T_prev.setZero();

        for (int iter = 0; iter < problem.getMaxIterations(); ++iter) {
            logger.info("--> Iteration ", iter + 1, " / ", problem.getMaxIterations());

            // --- Step 1: Solve EMag Field ---
            logger.info("    Solving EMag Field...");
            emag_field->assemble();
            emag_field->applyBCs();

            // Manually stabilize the matrix for the Temperature DOFs
            auto &emag_K = emag_field->getStiffnessMatrix();
            for (const auto &node: mesh.getNodes()) {
                int dof_idx = dof_manager.getEquationIndex(node->getId(), "Temperature");
                if (dof_idx != -1) emag_K.coeffRef(dof_idx, dof_idx) = 1.0;
            }

            LinearSolver::solve(emag_K, emag_field->getRHS(), emag_field->getSolution());

            // --- Step 2: Execute Coupling ---
            coupling_manager.executeCouplings();

            // --- Step 3: Solve Heat Field ---
            logger.info("    Solving Heat Field...");
            heat_field->assemble();
            heat_field->applyBCs();

            // Manually stabilize the matrix for the Voltage DOFs
            auto &heat_K = heat_field->getStiffnessMatrix();
            for (const auto &node: mesh.getNodes()) {
                int dof_idx = dof_manager.getEquationIndex(node->getId(), "Voltage");
                if (dof_idx != -1) heat_K.coeffRef(dof_idx, dof_idx) = 1.0;
            }

            LinearSolver::solve(heat_K, heat_field->getRHS(), heat_field->getSolution());

            // --- Step 4: Check for Convergence ---
            double norm_diff = (heat_field->getSolution() - T_prev).norm();
            double norm_sol = heat_field->getSolution().norm();
            double relative_error = (norm_sol > 1e-9) ? (norm_diff / norm_sol) : norm_diff;

            logger.info("    Convergence Check: Relative Error = ", relative_error);
            if (relative_error < problem.getConvergenceTolerance()) {
                logger.info("--- Coupled steady-state solver converged after ", iter + 1, " iterations. ---");
                return;
            }

            T_prev = heat_field->getSolution();
        }

        logger.warn("--- Coupled steady-state solver did not converge after ", problem.getMaxIterations(),
                    " iterations. ---");
    }

    void CoupledElectroThermalSolver::solveTransient(Core::Problem &problem) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("\n--- Starting Coupled Transient Solve ---");
        logger.info("Time Step: ", problem.getTimeStep(), "s, Total Time: ", problem.getTotalTime(), "s");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");
        auto &coupling_manager = problem.getCouplingManager();

        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * problem.getTimeStep(), "s");

            // 1. Solve EMag field based on previous temperature
            emag_field->assemble();
            emag_field->applyBCs();
            LinearSolver::solve(emag_field->getStiffnessMatrix(), emag_field->getRHS(), emag_field->getSolution());

            // 2. Execute coupling
            coupling_manager.executeCouplings();

            // 3. Solve Heat field for the current time step
            heat_field->assemble();
            Eigen::SparseMatrix<double> A = (heat_field->getMassMatrix() / problem.getTimeStep()) + heat_field->
                                            getStiffnessMatrix();
            Eigen::MatrixXd b = heat_field->getRHS() + (heat_field->getMassMatrix() / problem.getTimeStep()) *
                                heat_field->getPreviousSolution();

            auto A_bc = A;
            auto b_bc = b;
            for (const auto &bc: heat_field->getBCs()) {
                bc->apply(A_bc, b_bc);
            }
            LinearSolver::solve(A_bc, b_bc, heat_field->getSolution());

            // 4. Update temperature for the next step
            heat_field->updatePreviousSolution();
        }
        logger.info("\n--- Transient Solve Finished ---");
    }
} // Solver
