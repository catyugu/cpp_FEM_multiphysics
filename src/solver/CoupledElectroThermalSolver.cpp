//
// Created by HUAWEI on 2025/7/11.
//

#include "solver/CoupledElectroThermalSolver.hpp"

#include <physics/Current1D.hpp>
#include <physics/Heat1D.hpp>

#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "core/LinearSolver.hpp"
#include "utils/SimpleLogger.hpp"

using namespace Core;

namespace Solver {
    void CoupledElectroThermalSolver::solveSteadyState(Problem& problem) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("\n--- Solving Coupled Electro-Thermal Problem ---");

        auto *emag_field_base = dynamic_cast<Physics::Current2D *>(problem.getField("Voltage"));
        auto *heat_field_base = dynamic_cast<Physics::Heat2D *>(problem.getField("Temperature"));
        if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field_base)) {
            if (auto *heat1d = dynamic_cast<Physics::Heat1D *>(heat_field_base)) {
                emag1d->setCoupledHeatField(heat_field_base);
            }else {
                logger.error("Dimensions of the fields must be the same!");
                return;
            }
        } else if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field_base)) {
            if (auto *heat1d = dynamic_cast<Physics::Heat2D *>(heat_field_base)) {
                emag2d->setCoupledHeatField(heat_field_base);
            }else {
                logger.error("Dimensions of the fields must be the same!");
                return;
            }
        } else {
            logger.error("Unsupported field type for Coupled Electro-Thermal Solver.");
            return;
        }
        // 1. Ensure the heat field is unlinked before the first EMag solve.
        //    This forces EMag2D::assemble to use a default temperature, preventing a crash.
        emag_field_base->setCoupledHeatField(nullptr);

        // 2. Solve the EMag field first.
        logger.info("--> Solving EMag Field (assuming constant properties)...");
        emag_field_base->assemble();
        emag_field_base->applyBCs();
        LinearSolver::solve(emag_field_base->getStiffnessMatrix(), emag_field_base->getRHS(), emag_field_base->getSolution());

        // 3. Calculate Joule heat based on the solved voltage distribution
        //    and set it as the source for the thermal problem.
        auto joule_heat = emag_field_base->calculateJouleHeat();
        heat_field_base->setVolumetricHeatSource(joule_heat);

        // 4. Solve the Heat field with the calculated heat source.
        logger.info("--> Solving Heat Field...");
        heat_field_base->assemble();
        heat_field_base->applyBCs();
        LinearSolver::solve(heat_field_base->getStiffnessMatrix(), heat_field_base->getRHS(), heat_field_base->getSolution());
    }

    void CoupledElectroThermalSolver::solveTransient(Problem& problem) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("\n--- Starting Transient Solve ---");
        logger.info("Time Step: ", problem.getTimeStep(), "s, Total Time: ", problem.getTotalTime(), "s");

        auto *emag_field_base = problem.getField("Voltage");
        auto *heat_field_base = problem.getField("Temperature");

        if (emag_field_base && heat_field_base) {
            logger.info("\n--- Solving Coupled Transient Problem ---");

            // --- FIX: Implement the full coupled transient loop ---
            if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field_base)) {
                if (auto *heat1d = dynamic_cast<Physics::Heat1D *>(heat_field_base)) {
                    emag1d->setCoupledHeatField(heat_field_base);
                }else {
                    logger.error("Dimensions of the fields must be the same!");
                    return;
                }
            } else if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field_base)) {
                if (auto *heat1d = dynamic_cast<Physics::Heat2D *>(heat_field_base)) {
                    emag2d->setCoupledHeatField(heat_field_base);
                }else {
                    logger.error("Dimensions of the fields must be the same!");
                    return;
                }
            } else {
                logger.error("Unsupported field type for Coupled Electro-Thermal Solver.");
                return;
            }

            int num_steps = static_cast<int>(problem.getTotalTime() / problem.getTimeStep());
            for (int i = 0; i < num_steps; ++i) {
                logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i + 1) * problem.getTimeStep(), "s");

                emag_field_base->assemble();
                emag_field_base->applyBCs();
                LinearSolver::solve(emag_field_base->getStiffnessMatrix(), emag_field_base->getRHS(),
                                    emag_field_base->getSolution());

                std::vector<double> joule_heat;
                if (auto *emag1d = dynamic_cast<Physics::Current1D *>(emag_field_base)) {
                    joule_heat = emag1d->calculateJouleHeat();
                } else if (auto *emag2d = dynamic_cast<Physics::Current2D *>(emag_field_base)) {
                    joule_heat = emag2d->calculateJouleHeat();
                }

                if (auto *heat1d = dynamic_cast<Physics::Heat1D *>(heat_field_base)) {
                    heat1d->setVolumetricHeatSource(joule_heat);
                } else if (auto *heat2d = dynamic_cast<Physics::Heat2D *>(heat_field_base)) {
                    heat2d->setVolumetricHeatSource(joule_heat);
                }

                heat_field_base->assemble();
                Eigen::SparseMatrix<double> A = (heat_field_base->getMassMatrix() / problem.getTimeStep()) + heat_field_base->
                                                getStiffnessMatrix();
                Eigen::MatrixXd b = heat_field_base->getRHS() + (heat_field_base->getMassMatrix() / problem.getTimeStep()) *
                                    heat_field_base->getPreviousSolution();
                auto A_bc = A;
                auto b_bc = b;
                for (const auto &bc: heat_field_base->getBCs()) {
                    bc->apply(A_bc, b_bc);
                }
                LinearSolver::solve(A_bc, b_bc, heat_field_base->getSolution());
                heat_field_base->updatePreviousSolution();
            }

        }
    }


} // Solver