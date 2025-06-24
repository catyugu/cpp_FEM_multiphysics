#ifndef PHYSICSFIELD_HPP
#define PHYSICSFIELD_HPP

#include "core/Mesh.hpp"
#include "core/DOFManager.hpp"
#include "core/BoundaryCondition.hpp"
#include <vector>
#include <memory>
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace Physics {

    // Abstract base class representing a physical field (e.g., Heat, EM)
    class PhysicsField {
    public:
        virtual ~PhysicsField() = default;

        // Get the name of the physics field
        virtual const char* getName() const = 0;

        // Get the name of the variable this field solves for (e.g. "Voltage", "Temperature")
        virtual const char* getVariableName() const = 0;

        // Setup the field for a given mesh and DOFManager
        virtual void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) = 0;

        // Assemble the global matrices and vectors for the field
        virtual void assemble() = 0;

        // Add a boundary condition
        void addBC(std::unique_ptr<Core::BoundaryCondition> bc);

        // Apply all registered boundary conditions
        void applyBCs();

        // Getters for the system components
        Eigen::SparseMatrix<double>& getStiffnessMatrix();
        Eigen::VectorXd& getRHSVector();
        Eigen::VectorXd& getSolution();


    protected:
        Core::Mesh* mesh_ = nullptr;
        Core::DOFManager* dof_manager_ = nullptr;

        Eigen::SparseMatrix<double> K_; // Global stiffness matrix
        Eigen::VectorXd F_;               // Global force/RHS vector
        Eigen::VectorXd U_;               // Solution vector

        std::vector<std::unique_ptr<Core::BoundaryCondition>> bcs_;
    };

} // namespace Physics

#endif // PHYSICSFIELD_HPP
