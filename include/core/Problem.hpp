#ifndef PROBLEM_HPP
#define PROBLEM_HPP

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

        /**
         * @brief Sets the parameters for the iterative coupled solver.
         * @param max_iter The maximum number of iterations to perform.
         * @param tol The relative tolerance for checking convergence.
         */
        void setSolverParameters(int max_iter, double tol);

        void setup();
        void solve();
        void exportResults(const std::string& filename) const;

        // Const-correct getters
        Physics::PhysicsField* getField(const std::string& var_name) const;
        const Mesh& getMesh() const;
        const DOFManager& getDofManager() const;

    private:
        std::unique_ptr<Mesh> mesh_;
        std::unique_ptr<DOFManager> dof_manager_;
        std::vector<std::unique_ptr<Physics::PhysicsField>> fields_;

        // Solver control parameters
        int max_iterations_ = 20;
        double convergence_tolerance_ = 1e-4;
    };

} // namespace Core

#endif // PROBLEM_HPP
