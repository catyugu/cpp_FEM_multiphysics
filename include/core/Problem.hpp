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
        void setup();
        void solve();

        // Export all results to a file (e.g., VTK)
        void exportResults(const std::string& filename) const;

        // Const-correct getters
        Physics::PhysicsField* getField(const std::string& var_name) const;
        const Mesh& getMesh() const;
        const DOFManager& getDofManager() const;

    private:
        std::unique_ptr<Mesh> mesh_;
        std::unique_ptr<DOFManager> dof_manager_;
        std::vector<std::unique_ptr<Physics::PhysicsField>> fields_;
    };

} // namespace Core

#endif // PROBLEM_HPP
