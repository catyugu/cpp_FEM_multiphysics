//
// Created by HUAWEI on 2025/7/11.
//

#include <solver/SolverFactory.hpp>

#include "solver/SolverFactory.hpp"
#include "solver/SingleFieldSolver.hpp"
#include "solver/CoupledElectroThermalSolver.hpp"
#include "physics/PhysicsField.hpp"

namespace Solver {

    std::unique_ptr<Solver> SolverFactory::createSolver(const Core::Problem& problem) {
        auto *emag_field_base = problem.getField("Voltage");
        auto *heat_field_base = problem.getField("Temperature");


        if (emag_field_base && heat_field_base) {
            return std::make_unique<CoupledElectroThermalSolver>();
        } else {
            return std::make_unique<SingleFieldSolver>();
        }
    }

} // namespace Core
