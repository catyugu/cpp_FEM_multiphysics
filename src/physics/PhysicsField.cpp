#include "physics/PhysicsField.hpp"
#include "utils/SimpleLogger.hpp"

namespace Physics {

    void PhysicsField::addBC(std::unique_ptr<Core::BoundaryCondition> bc) {
        bcs_.push_back(std::move(bc));
    }

    void PhysicsField::applyBCs() {
        SimpleLogger::Logger::instance().info("Applying ", bcs_.size(), " BCs for ", getName());
        for(const auto& bc : bcs_) {
            bc->apply(K_, F_);
        }
    }

    // --- FIX: Implement new public accessors ---
    const std::vector<std::unique_ptr<Core::BoundaryCondition>>& PhysicsField::getBCs() const {
        return bcs_;
    }

    void PhysicsField::updatePreviousSolution() {
        U_prev_ = U_;
    }


    const Eigen::VectorXd& PhysicsField::getPreviousSolution() const {
        return U_prev_;
    }

    void PhysicsField::setInitialConditions(double initial_value) {
        if (U_.size() > 0) {
            U_.setConstant(initial_value);
            U_prev_ = U_;
            SimpleLogger::Logger::instance().info("Set initial condition for '", getVariableName(), "' to ", initial_value);
        } else {
            SimpleLogger::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(), "'.");
        }
    }

    Eigen::SparseMatrix<double>& PhysicsField::getStiffnessMatrix() { return K_; }
    Eigen::SparseMatrix<double>& PhysicsField::getMassMatrix() { return M_; }
    Eigen::VectorXd& PhysicsField::getRHSVector() { return F_; }
    Eigen::VectorXd& PhysicsField::getSolution() { return U_; }

    const Eigen::SparseMatrix<double>& PhysicsField::getStiffnessMatrix() const { return K_; }
    const Eigen::SparseMatrix<double>& PhysicsField::getMassMatrix() const { return M_; }
    const Eigen::VectorXd& PhysicsField::getRHSVector() const { return F_; }
    const Eigen::VectorXd& PhysicsField::getSolution() const { return U_; }

} // namespace Physics
