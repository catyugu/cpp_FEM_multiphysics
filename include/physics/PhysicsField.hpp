#ifndef PHYSICSFIELD_HPP
#define PHYSICSFIELD_HPP

#include <core/mesh/Mesh.hpp>
#include "core/DOFManager.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include <vector>
#include <memory>
#include <core/material/Material.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <core/sources/SourceTerm.hpp>
#include <functional>

// Forward declarations
namespace Core {
    class Problem;
}

namespace Post { class HeatFluxCalculator; }
namespace IO { class Exporter; }
namespace Solver { class CoupledElectroThermalSolver; }
namespace Core { class ElectroThermalCoupling; }

namespace Physics {
    class PhysicsField {
        friend class Core::ElectroThermalCoupling;
        friend class Solver::CoupledElectroThermalSolver;
        friend class IO::Exporter;
        friend class Post::HeatFluxCalculator;

    public:
        virtual ~PhysicsField() = default;

        virtual const char *getName() const = 0;
        virtual const char *getVariableName() const = 0;
        virtual int getDimension() const = 0;
        virtual int getNumComponents() const { return 1; }
        virtual const Core::Material& getMaterial(const Core::Element* elem) const = 0;

        virtual void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager);
        virtual void assemble() = 0;

        void addBC(std::unique_ptr<Core::BoundaryCondition> bc);
        void addBCs(std::vector<std::unique_ptr<Core::BoundaryCondition>> &&bcs);
        void applyBCs();
        const std::vector<std::unique_ptr<Core::BoundaryCondition> > &getBCs() const;

        void addSource(std::unique_ptr<Core::SourceTerm> source);
        void applySources();

        void updatePreviousSolution();
        const Eigen::VectorXd &getPreviousSolution() const;
        void setInitialConditions(double initial_value);
        template<typename Func>
        void setInitialConditions(Func func);

        void setElementOrder(int order);
        int getElementOrder() const { return element_order_; }

        Eigen::SparseMatrix<double> &getStiffnessMatrix();
        Eigen::SparseMatrix<double> &getMassMatrix();
        Eigen::VectorXd &getRHS();
        Eigen::VectorXd &getCouplingRHS(); // 获取耦合源向量
        Eigen::VectorXd &getSolution();

        const Eigen::SparseMatrix<double> &getStiffnessMatrix() const;
        const Eigen::SparseMatrix<double> &getMassMatrix() const;
        const Eigen::VectorXd &getRHS() const;
        const Eigen::VectorXd &getSolution() const;

        const Core::Mesh *getMesh() const { return mesh_; }
        const Core::DOFManager *getDofManager() const { return dof_manager_; }
        std::vector<int> getElementDofs(Core::Element *elem) const;

        void removeBCsByTag(const std::string &tag);
        void removeSourcesByTag(const std::string &tag);

        // 新增获取方法
        const std::vector<std::unique_ptr<Core::SourceTerm>>& getSources() const { return source_terms_; }

        void enable() { enabled = true; }
        void disable() { enabled = false; }
        bool isEnabled() const { return enabled; }

        // ========= 新增：变量管理和更新功能 =========

        /**
         * @brief Update element variable values after solution
         * This method should be called after solving to update all element variable values
         * based on the current solution vector
         */
        virtual void updateElementVariables();

        /**
         * @brief Get the variable name associated with this physics field
         * This allows automatic variable registration and updates
         * @return The variable name (e.g., "T" for temperature, "V" for voltage)
         */
        virtual const std::string& getFieldVariableName() const = 0;

        /**
         * @brief Register this field's variable with the global variable manager
         * This should be called during setup to ensure the variable exists
         */
        void registerFieldVariable();

    protected:
        Core::Problem* problem_ = nullptr;
        Core::Mesh *mesh_ = nullptr;
        Core::DOFManager *dof_manager_ = nullptr;
        bool enabled = true;
        int element_order_ = 1;

        Eigen::SparseMatrix<double> K_;
        Eigen::SparseMatrix<double> M_;
        Eigen::VectorXd F_;
        Eigen::VectorXd F_coupling_;
        Eigen::VectorXd U_;
        Eigen::VectorXd U_prev_;

        std::vector<std::unique_ptr<Core::BoundaryCondition> > bcs_;
        std::vector<std::unique_ptr<Core::SourceTerm> > source_terms_;
    };

}

#endif // PHYSICSFIELD_HPP

