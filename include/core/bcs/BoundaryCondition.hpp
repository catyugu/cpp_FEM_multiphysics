#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <string>
#include <core/mesh/Mesh.hpp>

// Forward declaration
namespace Core {
    class DOFManager;
}

namespace Core {

    class BoundaryCondition {
    public:
        explicit BoundaryCondition(std::string tag = "") : tag_(std::move(tag)) {}
        virtual ~BoundaryCondition() = default;

        virtual void apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const = 0;
        virtual int getEquationIndex() const = 0;
        const std::string& getTag() const { return tag_; }
    protected:
        std::string tag_;
    };

    class DirichletBC : public BoundaryCondition {
    public:
        DirichletBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd value, const std::string& tag = "");
        DirichletBC(int equation_index, Eigen::VectorXd value, const std::string& tag = "");
        static std::vector<std::unique_ptr<BoundaryCondition>> create(
            const DOFManager& dof_manager,
            const Mesh& mesh,
            const std::string& var_name,
            int field_order,
            const std::function<bool(const std::vector<double>& coords)>& region_predicate,
            const std::function<double(const std::vector<double>& coords)>& value_function,
            const std::string& tag = ""
        );

        void apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const override;
        int getEquationIndex() const override { return equation_index_; }
        const Eigen::VectorXd& getValue() const { return value_; }

    private:
        int equation_index_;
        Eigen::VectorXd value_;
    };


    class NeumannBC : public BoundaryCondition {
    public:
        NeumannBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd flux_value, const std::string& tag = "");
        void apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const override;
        int getEquationIndex() const override { return equation_index_; }
    private:
        int equation_index_;
        Eigen::VectorXd flux_value_;
    };

    class CauchyBC : public BoundaryCondition {
    public:
        CauchyBC(const DOFManager& dof_manager, int node_id, const std::string& var_name, Eigen::VectorXd h, Eigen::VectorXd T_inf, const std::string& tag = "");
        void apply(Eigen::SparseMatrix<double>& K, Eigen::VectorXd& F) const override;
        int getEquationIndex() const override { return equation_index_; }
    private:
        int equation_index_;
        Eigen::VectorXd h_;
        Eigen::VectorXd T_inf_;
    };

} // namespace Core

#endif // BOUNDARYCONDITION_HPP