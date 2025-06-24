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

    class PhysicsField {
    public:
        virtual ~PhysicsField() = default;

        virtual const char* getName() const = 0;
        virtual const char* getVariableName() const = 0;

        virtual void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) = 0;
        virtual void assemble() = 0;

        void addBC(std::unique_ptr<Core::BoundaryCondition> bc);
        void applyBCs();

        // --- FIX: Add public accessors for transient solver ---
        const std::vector<std::unique_ptr<Core::BoundaryCondition>>& getBCs() const;
        void updatePreviousSolution();

        // Getter for the solution at the previous time step
        const Eigen::VectorXd& getPreviousSolution() const;
        void setInitialConditions(double initial_value);

        // Getters for matrices and vectors
        Eigen::SparseMatrix<double>& getStiffnessMatrix();
        Eigen::SparseMatrix<double>& getMassMatrix();
        Eigen::VectorXd& getRHSVector();
        Eigen::VectorXd& getSolution();

        const Eigen::SparseMatrix<double>& getStiffnessMatrix() const;
        const Eigen::SparseMatrix<double>& getMassMatrix() const;
        const Eigen::VectorXd& getRHSVector() const;
        const Eigen::VectorXd& getSolution() const;

    protected:
        Core::Mesh* mesh_ = nullptr;
        Core::DOFManager* dof_manager_ = nullptr;

        Eigen::SparseMatrix<double> K_;
        Eigen::SparseMatrix<double> M_;
        Eigen::VectorXd F_;
        Eigen::VectorXd U_;
        Eigen::VectorXd U_prev_;

        std::vector<std::unique_ptr<Core::BoundaryCondition>> bcs_;
    };

} // namespace Physics

#endif // PHYSICSFIELD_HPP
