//
// Created by HUAWEI on 2025/7/11.
//

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "core/Problem.hpp"

namespace Solver {
    class Solver {
    public:
        virtual ~Solver() = default;
        virtual void solveSteadyState(Core::Problem& problem) = 0;
        virtual void solveTransient(Core::Problem& problem) = 0;
    };
}

#endif //SOLVER_HPP
