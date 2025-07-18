#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <map>
#include <vector>
#include <memory>
#include <string>

#include <core/coupling/CouplingManager.hpp>
#include <solver/Solver.hpp>
#include <solver/LinearSolver.hpp>
#include "post/PostProcessor.hpp" // Include the new PostProcessor header

// Forward declarations
namespace Core {
    class Mesh;
    class DOFManager;
}
namespace Physics {
    class PhysicsField;
}
namespace Solver {
    class Solver;
}

namespace Core {

    class Problem {
    public:
        explicit Problem(std::unique_ptr<Mesh> mesh);
        ~Problem();

        void addField(std::unique_ptr<Physics::PhysicsField> field);
        void addPostProcessor(std::unique_ptr<Post::PostProcessor> post_processor); // New method
        void setup();

        // --- Solver Control ---
        void setIterativeSolverParameters(int max_iter, double tol);
        void setTimeSteppingParameters(double time_step, double total_time);
        void setLinearSolverType(Solver::SolverType type);

        // --- Solution Methods ---
        void solveSteadyState();
        void solveTransient();

        void exportResults(const std::string& filename) const;
        const Mesh& getMesh() const;
        const DOFManager& getDofManager() const;
        Physics::PhysicsField* getField(const std::string& var_name) const;
        const std::vector<std::unique_ptr<Physics::PhysicsField>>& getFields() const;
        const std::map<std::string, Post::PostProcessingResult>& getPostProcessingResults() const; // New accessor

        // Solver control parameters
        int getMaxIterations() const { return max_iterations_;}
        double getConvergenceTolerance() const { return convergence_tolerance_;}
        double getTimeStep() const { return time_step_;}
        double getTotalTime() const { return total_time_;}
        Solver::SolverType getLinearSolverType() const { return linear_solver_type_; }

        CouplingManager& getCouplingManager() { return coupling_manager_; }

    private:
        void runPostProcessors(); // New private helper method

        std::unique_ptr<Mesh> mesh_;
        std::unique_ptr<DOFManager> dof_manager_;
        std::unique_ptr<Solver::Solver> solver_;
        std::vector<std::unique_ptr<Physics::PhysicsField>> fields_;
        std::vector<std::unique_ptr<Post::PostProcessor>> post_processors_; // New member
        std::map<std::string, Post::PostProcessingResult> post_processing_results_; // New member
        CouplingManager coupling_manager_;

        int max_iterations_ = 20;
        double convergence_tolerance_ = 1e-4;
        double time_step_ = 0.1;
        double total_time_ = 1.0;
        Solver::SolverType linear_solver_type_ = Solver::SolverType::LU;
    };

} // namespace Core

#endif // PROBLEM_HPP
