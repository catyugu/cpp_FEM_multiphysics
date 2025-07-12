#include "solver/CoupledElectroThermalSolver.hpp"
#include "core/Problem.hpp"
#include "physics/Current1D.hpp"
#include "physics/Heat1D.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include <solver/LinearSolver.hpp>
#include "utils/SimpleLogger.hpp"

namespace Solver {
    void CoupledElectroThermalSolver::solveSteadyState(Core::Problem& problem) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem ---");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");

        // Uncouple for the first solve to use default temperature
        if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field)) {
            emag2d->setCoupledHeatField(nullptr);
        } else if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field)) {
            emag1d->setCoupledHeatField(nullptr);
        }

        logger.info("--> Solving EMag Field (initial)...");
        emag_field->assemble();
        emag_field->applyBCs();
        LinearSolver::solve(emag_field->getStiffnessMatrix(), emag_field->getRHS(), emag_field->getSolution());

        // Couple the fields for subsequent calculations
        if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field)) {
            emag2d->setCoupledHeatField(heat_field);
        } else if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field)) {
            emag1d->setCoupledHeatField(heat_field);
        }

        std::vector<double> joule_heat;
        if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field)) {
            joule_heat = emag1d->calculateJouleHeat();
            dynamic_cast<Physics::Heat1D*>(heat_field)->setVolumetricHeatSource(joule_heat);
        } else if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field)) {
            joule_heat = emag2d->calculateJouleHeat();
            dynamic_cast<Physics::Heat2D*>(heat_field)->setVolumetricHeatSource(joule_heat);
        }

        logger.info("--> Solving Heat Field...");
        heat_field->assemble();
        heat_field->applyBCs();
        LinearSolver::solve(heat_field->getStiffnessMatrix(), heat_field->getRHS(), heat_field->getSolution());
    }

    void CoupledElectroThermalSolver::solveTransient(Core::Problem& problem) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("\n--- Starting Coupled Transient Solve ---");
        logger.info("Time Step: ", problem.getTimeStep(), "s, Total Time: ", problem.getTotalTime(), "s");

        auto *emag_field = problem.getField("Voltage");
        auto *heat_field = problem.getField("Temperature");

        int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * problem.getTimeStep(), "s");

            // 1. Solve EMag field based on previous temperature
            emag_field->assemble();
            emag_field->applyBCs();
            LinearSolver::solve(emag_field->getStiffnessMatrix(), emag_field->getRHS(), emag_field->getSolution());

            // 2. Calculate Joule Heat and set as source for thermal problem
            std::vector<double> joule_heat;
            if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field)) {
                joule_heat = emag1d->calculateJouleHeat();
                dynamic_cast<Physics::Heat1D*>(heat_field)->setVolumetricHeatSource(joule_heat);
            } else if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field)) {
                joule_heat = emag2d->calculateJouleHeat();
                dynamic_cast<Physics::Heat2D*>(heat_field)->setVolumetricHeatSource(joule_heat);
            }

            // 3. Solve Heat field for the current time step
            heat_field->assemble();
            Eigen::SparseMatrix<double> A = (heat_field->getMassMatrix() / problem.getTimeStep()) + heat_field->getStiffnessMatrix();
            Eigen::MatrixXd b = heat_field->getRHS() + (heat_field->getMassMatrix() / problem.getTimeStep()) * heat_field->getPreviousSolution();

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