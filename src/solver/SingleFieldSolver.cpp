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
            field->assemble();
            field->applyBCs();
            LinearSolver::solve(field->getStiffnessMatrix(), field->getRHS(), field->getSolution());
        }
    }

    void SingleFieldSolver::solveTransient(Core::Problem &problem) {
        auto &logger = Utils::Logger::instance();
        logger.info("\n--- Solving Single-Field Transient Problem ---");
        auto* field = problem.getFields()[0].get(); // Assuming one field
        field->assemble();

        Eigen::SparseMatrix<double> A_eff;
        Eigen::MatrixXd b_eff;

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

            LinearSolver::solve(A_bc, b_bc, field->getSolution());
            field->updatePreviousSolution();
        }
    }
} // Solver