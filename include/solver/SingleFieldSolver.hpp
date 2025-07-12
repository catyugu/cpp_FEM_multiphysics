#ifndef SINGLEFIELDSOLVER_HPP
#define SINGLEFIELDSOLVER_HPP
#include "Solver.hpp"
namespace Solver {

    class SingleFieldSolver : public Solver {
        void solveTransient(Core::Problem &problem) override;
        void solveSteadyState(Core::Problem &problem) override;
    };

} // Solver

#endif //SINGLEFIELDSOLVER_HPP