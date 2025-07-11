//
// Created by HUAWEI on 2025/7/11.
//

#ifndef SINGLEFIELDSOLVER_HPP
#define SINGLEFIELDSOLVER_HPP
#include "Solver.hpp"
namespace Solver {

class SingleFieldSolver : public Solver {
    void solveSteadyState(Core::Problem& problem) override;
    void solveTransient(Core::Problem& problem) override;
};

} // Solver

#endif //SINGLEFIELDSOLVER_HPP
