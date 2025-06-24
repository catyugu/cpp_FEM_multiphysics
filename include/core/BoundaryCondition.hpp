#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include "DOFManager.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace Core {

    // Abstract base class for Boundary Conditions
    class BoundaryCondition {
    public:
        virtual ~BoundaryCondition() = default;

        // Apply the boundary condition to the global system of equations
        virtual void apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const = 0;
    };

    // Dirichlet (fixed value) Boundary Condition
    class DirichletBC : public BoundaryCondition {
    public:
        DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double value);

        void apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const override;

    private:
        int equation_index_;
        double value_;
    };

} // namespace Core

#endif // BOUNDARYCONDITION_HPP
