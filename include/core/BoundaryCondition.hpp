#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <string>
#include <memory>

// Forward declaration to prevent circular includes
namespace Core {
    class DOFManager;
}

namespace Core {

// Abstract base class for all Boundary Conditions
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    // Apply the boundary condition to the global system of equations
    virtual void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const = 0;
};

// --- Concrete BC Implementations ---

// Dirichlet (Type 1): Specifies a fixed value for a DOF (e.g., T = 373K)
class DirichletBC : public BoundaryCondition {
public:
    DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double value);
    void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override;

private:
    int equation_index_;
    double value_;
};

// Neumann (Type 2): Specifies a flux at a node (e.g., heat flux q = 100 W/m^2)
// For 1D, this is a point value. In 2D/3D, this would be integrated over an edge/face.
class NeumannBC : public BoundaryCondition {
public:
    NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double flux_value);
    void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override;

private:
    int equation_index_;
    double flux_value_;
};

// Cauchy/Robin/Mixed (Type 3): Specifies a relationship between a value and its flux (e.g., convection)
// h * (T_inf - T_surface) = -k * dT/dn
// This adds to both the stiffness matrix K and the RHS vector F.
class CauchyBC : public BoundaryCondition {
public:
    CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, double h, double T_inf);
    void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override;

private:
    int equation_index_;
    double h_;     // Convection coefficient
    double T_inf_; // Ambient temperature
};


} // namespace Core

#endif // BOUNDARYCONDITION_HPP
