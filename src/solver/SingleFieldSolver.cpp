#include "solver/SingleFieldSolver.hpp"
#include "core/Problem.hpp"
#include <solver/LinearSolver.hpp>
#include "physics/PhysicsField.hpp"
#include "utils/SimpleLogger.hpp"

namespace Solver {
    void SingleFieldSolver::solveSteadyState(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Uncoupled Steady-State Physics ---");
        for (const auto &field: problem.getFields()) {
            if (!field->isEnabled()) {
                continue;
            }
            logger.info("Solving for field: ", field->getName());
            field->assemble(nullptr);
            field->applySources();
            field->applyBCs();

            LinearSolver::solve(field->getStiffnessMatrix(), field->getRHS(), field->getSolution(),
                                problem.getLinearSolverType(),
                                problem.getMaxIterations(),
                                problem.getConvergenceTolerance());
        }
    }

    void SingleFieldSolver::solveTransient(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Single-Field Transient Problem ---");

        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        double dt = problem.getTimeStep();

        for (const auto& field : problem.getFields()) {
            if (!field->isEnabled()) continue;
            field->updatePreviousSolution(); // Initialize U_prev_ from U_
        }

        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * dt, "s");

            for (const auto& field : problem.getFields()) {
                if (!field->isEnabled()) continue;

                logger.info("  -> Assembling for field: ", field->getName());
                field->assemble(nullptr);
                field->applySources(); // Important to apply time-dependent sources inside the loop

                // Backward Euler scheme: (M/dt + K) * U_n+1 = F + (M/dt) * U_n
                Eigen::SparseMatrix<double> A_eff = (field->getMassMatrix() / dt) + field->getStiffnessMatrix();
                Eigen::VectorXd b_eff = field->getRHS() + (field->getMassMatrix() / dt) * field->getPreviousSolution();

                // Apply BCs to the effective system
                for (const auto &bc: field->getBCs()) {
                    bc->apply(A_eff, b_eff);
                }

                LinearSolver::solve(A_eff, b_eff, field->getSolution(),
                                    problem.getLinearSolverType(),
                                    problem.getMaxIterations(),
                                    problem.getConvergenceTolerance());

                field->updatePreviousSolution(); // U_prev_ becomes the solution we just found
            }
        }
    }
} // Solver