//
// Created by HUAWEI on 2025/7/11.
//

#ifndef SOLVERFACTORY_HPP
#define SOLVERFACTORY_HPP
#include <memory>
#include <core/Problem.hpp>
#include <solver/Solver.hpp>

namespace Solver{

    class SolverFactory {
    public:
        static std::unique_ptr<Solver> createSolver(const Core::Problem& problem);
    };

}



#endif //SOLVERFACTORY_HPP
