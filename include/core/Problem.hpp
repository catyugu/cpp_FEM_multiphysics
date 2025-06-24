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

    // Orchestrates the setup and solution of a multiphysics problem.
    class Problem {
    public:
        explicit Problem(std::unique_ptr<Mesh> mesh);
        ~Problem(); // Add destructor for unique_ptr forward decl

        // Add a physics field to the problem
        void addField(std::unique_ptr<Physics::PhysicsField> field);

        // Initialize all components
        void setup();

        // Run the solver workflow
        void solve();

        // Get a pointer to a specific field by its variable name
        Physics::PhysicsField* getField(const std::string& var_name);

        // Getters for core components
        Mesh& getMesh() const;
        DOFManager& getDofManager() const;

    private:
        std::unique_ptr<Mesh> mesh_;
        std::unique_ptr<DOFManager> dof_manager_;
        std::vector<std::unique_ptr<Physics::PhysicsField>> fields_;
    };

} // namespace Core

#endif // PROBLEM_HPP
