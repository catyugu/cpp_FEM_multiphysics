#ifndef SOLVER_HPP
#define SOLVER_HPP

namespace Core {
    class Problem;
}

namespace Solver {
    class Solver {
    public:
        virtual ~Solver() = default;
        virtual void solveSteadyState(Core::Problem& problem) = 0;
        virtual void solveTransient(Core::Problem& problem) = 0;
    };
}

#endif //SOLVER_HPP