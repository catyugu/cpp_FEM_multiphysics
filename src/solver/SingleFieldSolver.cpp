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
            field->applySources(); // Apply sources AFTER assembly
            field->applyBCs();

            LinearSolver::solve(field->getStiffnessMatrix(), field->getRHS(), field->getSolution(),
                                problem.getLinearSolverType(),
                                problem.getMaxIterations(),
                                problem.getConvergenceTolerance());
        }
    }

    // The solveTransient method also needs to be updated for consistency
    void SingleFieldSolver::solveTransient(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Single-Field Transient Problem ---");
        auto* field = problem.getFields()[0].get(); // Assuming one field
        field->assemble(nullptr);

        Eigen::SparseMatrix<double> A_eff;
        Eigen::VectorXd b_eff;

        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * problem.getTimeStep(), "s");
            A_eff = (field->getMassMatrix() / problem.getTimeStep()) + field->getStiffnessMatrix();
            b_eff = field->getRHS() + (field->getMassMatrix() / problem.getTimeStep()) * field->getPreviousSolution();

            auto A_bc = A_eff;
            auto b_bc = b_eff;
            for (const auto &bc: field->getBCs()) {
                bc->apply(A_bc, b_bc);
            }

            // Also apply the fix here
            LinearSolver::solve(A_bc, b_bc, field->getSolution(),
                                problem.getLinearSolverType(),
                                problem.getMaxIterations(),
                                problem.getConvergenceTolerance());

            field->updatePreviousSolution();
        }
    }
} // Solver