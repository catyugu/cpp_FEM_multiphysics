#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <string>

// Forward declaration
namespace Core {
    class DOFManager;
}

namespace Core {

    // Abstract base class for all Boundary Conditions
    class BoundaryCondition {
    public:
        explicit BoundaryCondition(std::string tag = "") : tag_(std::move(tag)) {}
        virtual ~BoundaryCondition() = default;

        virtual void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const = 0;
        const std::string& getTag() const { return tag_; }

    protected:
        std::string tag_;
    };

    // --- Concrete BC Implementations ---

    class DirichletBC : public BoundaryCondition {
    public:
        DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value, const std::string& tag = "");
        void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override;

    private:
        int equation_index_;
        Eigen::VectorXd value_;
    };

    class NeumannBC : public BoundaryCondition {
    public:
        NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value, const std::string& tag = "");
        void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override;

    private:
        int equation_index_;
        Eigen::VectorXd flux_value_;
    };

    class CauchyBC : public BoundaryCondition {
    public:
        CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf, const std::string& tag = "");
        void apply(Eigen::SparseMatrix<double>& K, Eigen::MatrixXd& F) const override;

    private:
        int equation_index_;
        Eigen::VectorXd h_;
        Eigen::VectorXd T_inf_;
    };

} // namespace Core

#endif // BOUNDARYCONDITION_HPP