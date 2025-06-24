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

    // --- Non-const Getters ---
    Eigen::SparseMatrix<double>& PhysicsField::getStiffnessMatrix() { return K_; }
    Eigen::VectorXd& PhysicsField::getRHSVector() { return F_; }
    Eigen::VectorXd& PhysicsField::getSolution() { return U_; }

    // --- FIX: Add implementations for const Getters ---
    const Eigen::SparseMatrix<double>& PhysicsField::getStiffnessMatrix() const { return K_; }
    const Eigen::VectorXd& PhysicsField::getRHSVector() const { return F_; }
    const Eigen::VectorXd& PhysicsField::getSolution() const { return U_; }

} // namespace Physics
