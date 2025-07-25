#ifndef PHYSICSFIELD_HPP
#define PHYSICSFIELD_HPP

#include <core/mesh/Mesh.hpp>
#include "core/DOFManager.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include <vector>
#include <memory>
#include <core/Material.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <core/sources/SourceTerm.hpp>
#include <functional> // For std::function

namespace Post {
    class HeatFluxCalculator;
}

namespace IO {
    class Exporter;
}

namespace Solver {
    class CoupledElectroThermalSolver;
}

namespace Core {
    class ElectroThermalCoupling;
}

namespace Physics {

    class PhysicsField {
        friend class Core::ElectroThermalCoupling;
        friend class Solver::CoupledElectroThermalSolver;
        friend class IO::Exporter;
        friend class Post::HeatFluxCalculator;
    public:
        virtual ~PhysicsField() = default;

        virtual const char* getName() const = 0;
        virtual const char* getVariableName() const = 0;
        virtual const Core::Material& getMaterial() const = 0;

        virtual int getDimension() const = 0; // New virtual function for dimension

        virtual void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) = 0;
        virtual void assemble() = 0;

        void addBC(std::unique_ptr<Core::BoundaryCondition> bc);
        void addBCs(std::vector<std::unique_ptr<Core::BoundaryCondition>> bcs);
        void applyBCs();


        const std::vector<std::unique_ptr<Core::BoundaryCondition>>& getBCs() const;
        void updatePreviousSolution();

        const Eigen::VectorXd& getPreviousSolution() const;

        void setInitialConditions(double initial_value);
        template<typename F>
        void setInitialConditions(std::function<F> initial_conditions);

        void setElementOrder(int order);

        Eigen::SparseMatrix<double>& getStiffnessMatrix();
        Eigen::SparseMatrix<double>& getMassMatrix();
        Eigen::VectorXd& getRHS();
        Eigen::VectorXd& getSolution();


        const Eigen::SparseMatrix<double>& getStiffnessMatrix() const;
        const Eigen::SparseMatrix<double>& getMassMatrix() const;
        const Eigen::VectorXd& getRHS() const;
        const Eigen::VectorXd& getSolution() const;
        const Core::Mesh* getMesh() const { return mesh_;}
        // get DOFManager

        const Core::DOFManager* getDofManager() const { return dof_manager_;}
        void removeBCsByTag(const std::string& tag);
        int getElementOrder() const { return element_order_; }

        void addSource(std::unique_ptr<Core::SourceTerm> source);
        void removeSourcesByTag(const std::string& tag);
        void applySources();

        void enable(){ enabled = true;}
        void disable(){ enabled = false;}
        bool isEnabled() const { return enabled; }
    protected:
        std::vector<int> get_element_dofs(Core::Element* elem) const; // New function for getting DOFs
        Core::Mesh* mesh_ = nullptr;
        Core::DOFManager* dof_manager_ = nullptr;
        bool enabled = true;
        int element_order_ = 1;

        Eigen::SparseMatrix<double> K_;
        Eigen::SparseMatrix<double> M_;
        Eigen::VectorXd F_;
        Eigen::VectorXd U_;
        Eigen::VectorXd U_prev_;

        std::vector<std::unique_ptr<Core::BoundaryCondition>> bcs_;
        std::vector<std::unique_ptr<Core::SourceTerm>> source_terms_;
    };

} // namespace Physics

#endif // PHYSICSFIELD_HPP