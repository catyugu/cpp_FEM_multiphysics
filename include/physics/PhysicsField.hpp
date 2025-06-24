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

        // --- FIX: Provide both const and non-const overloads for getters ---
        // Non-const version for modification
        Eigen::SparseMatrix<double>& getStiffnessMatrix();
        Eigen::VectorXd& getRHSVector();
        Eigen::VectorXd& getSolution();

        // Const version for read-only access
        const Eigen::SparseMatrix<double>& getStiffnessMatrix() const;
        const Eigen::VectorXd& getRHSVector() const;
        const Eigen::VectorXd& getSolution() const;

    protected:
        Core::Mesh* mesh_ = nullptr;
        Core::DOFManager* dof_manager_ = nullptr;

        Eigen::SparseMatrix<double> K_;
        Eigen::VectorXd F_;
        Eigen::VectorXd U_;

        std::vector<std::unique_ptr<Core::BoundaryCondition>> bcs_;
    };

} // namespace Physics

#endif // PHYSICSFIELD_HPP
