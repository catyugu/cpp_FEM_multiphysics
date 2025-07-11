//
// Created by HUAWEI on 2025/7/11.
//

#ifndef COUPLEDELECTROTHERMALSOLVER_HPP
#define COUPLEDELECTROTHERMALSOLVER_HPP

#include <core/Problem.hpp>
#include <solver/Solver.hpp>

namespace Solver {

    class CoupledElectroThermalSolver : public Solver {
    public:
        void solveSteadyState(Core::Problem& problem) override;
        void solveTransient(Core::Problem& problem) override;
    };
} // Solver

#endif //COUPLEDELECTROTHERMALSOLVER_HPP
