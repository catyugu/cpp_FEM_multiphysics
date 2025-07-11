#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <map>
#include <vector>
#include <memory>
#include <string>

// Forward declarations
namespace Core {
    class Mesh;
    class DOFManager;
}
namespace Physics {
    class PhysicsField;
}

namespace Core {

    class Problem {
    public:
        explicit Problem(std::unique_ptr<Mesh> mesh);
        ~Problem();

        void addField(std::unique_ptr<Physics::PhysicsField> field);
        void setup();

        // --- Solver Control ---
        void setIterativeSolverParameters(int max_iter, double tol);
        void setTimeSteppingParameters(double time_step, double total_time);

        // --- Solution Methods ---
        void solveSteadyState() const;
        void solveTransient() const;

        void exportResults(const std::string& filename) const;
        const Mesh& getMesh() const;
        const DOFManager& getDofManager() const;
        Physics::PhysicsField* getField(const std::string& var_name) const;

        // get all the  fields
        const std::vector<std::unique_ptr<Physics::PhysicsField>>& getFields() const;

        // Solver control parameters
        int getMaxIterations() const { return max_iterations_;}
        double getConvergenceTolerance() const { return convergence_tolerance_;}
        double getTimeStep() const { return time_step_;}
        double getTotalTime() const { return total_time_;}

    private:
        std::unique_ptr<Mesh> mesh_;
        std::unique_ptr<DOFManager> dof_manager_;
        std::vector<std::unique_ptr<Physics::PhysicsField>> fields_;

        int max_iterations_ = 20;
        double convergence_tolerance_ = 1e-4;
        double time_step_ = 0.1;
        double total_time_ = 1.0;
    };

} // namespace Core

#endif // PROBLEM_HPP
