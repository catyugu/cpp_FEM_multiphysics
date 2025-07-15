#include "physics/PhysicsField.hpp"
#include "utils/SimpleLogger.hpp"

namespace Physics {

    void PhysicsField::addBC(std::unique_ptr<Core::BoundaryCondition> bc) {
        bcs_.push_back(std::move(bc));
    }

    void PhysicsField::applyBCs() {
        auto& logger = Utils::Logger::instance();
        logger.info("Applying ", bcs_.size(), " defined BCs for ", getName());

        std::map<int, const Core::BoundaryCondition*> consolidated_bcs;
        for (const auto& bc : bcs_) {
            if(bc->getEquationIndex() != -1) {
                consolidated_bcs[bc->getEquationIndex()] = bc.get();
            }
        }

        logger.info("Applying ", consolidated_bcs.size(), " unique, consolidated BCs.");

        // Iterate through the consolidated map and apply the final, unique BCs.
        for (const auto& pair : consolidated_bcs) {
            pair.second->apply(K_, F_);
        }
        // --- END OF FIX ---
    }


    void PhysicsField::applySources() {
        // This is the new implementation
        F_.setZero(); // Start with a fresh force vector
        for (const auto& source : source_terms_) {
            source->apply(F_, *dof_manager_, *mesh_, getVariableName());
        }
    }


    void PhysicsField::removeBCsByTag(const std::string& tag) {
        if (tag.empty()) return; // Do not remove BCs with no tag
        auto it = std::remove_if(bcs_.begin(), bcs_.end(),
            [&](const std::unique_ptr<Core::BoundaryCondition>& bc_ptr) {
                return bc_ptr->getTag() == tag;
            });

        if (it != bcs_.end()) {
            bcs_.erase(it, bcs_.end());
        }
    }
    void PhysicsField::addSource(std::unique_ptr<Core::SourceTerm> source) {
        source_terms_.push_back(std::move(source));
    }

    void PhysicsField::removeSourcesByTag(const std::string& tag) {
        auto it = std::remove_if(source_terms_.begin(), source_terms_.end(),
            [&](const std::unique_ptr<Core::SourceTerm>& source_ptr) {
                return source_ptr->getTag() == tag;
            });

        if (it != source_terms_.end()) {
            source_terms_.erase(it, source_terms_.end());
        }
    }
    const std::vector<std::unique_ptr<Core::BoundaryCondition>>& PhysicsField::getBCs() const {
        return bcs_;
    }

    void PhysicsField::updatePreviousSolution() {
        U_prev_ = U_;
    }

    void PhysicsField::setElementOrder(int order) {
        if (order < 1) {
            Utils::Logger::instance().error("Element order must be at least 1.");
            return;
        }
        element_order_ = order;
    }


    const Eigen::MatrixXd& PhysicsField::getPreviousSolution() const {
        return U_prev_;
    }

    void PhysicsField::setInitialConditions(double initial_value) {
        if (U_.size() > 0) {
            U_.setConstant(initial_value);
            U_prev_ = U_;
            Utils::Logger::instance().info("Set initial condition for '", getVariableName(), "' to ", initial_value);
        } else {
            Utils::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(), "'.");
        }
    }

    template<typename F>
    void PhysicsField::setInitialConditions(std::function<F> func) {
        if (U_.size() > 0) {
            for (int i = 0; i < U_.size(); ++i) {
                U_(i) = func(mesh_->getNode(i));
            }
            U_prev_ = U_;
            Utils::Logger::instance().info("Set initial condition for '", getVariableName(), "' to ", func);
        } else {
            Utils::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(), "'.");
        }
    }


    Eigen::SparseMatrix<double>& PhysicsField::getStiffnessMatrix() { return K_; }
    Eigen::SparseMatrix<double>& PhysicsField::getMassMatrix() { return M_; }
    Eigen::MatrixXd& PhysicsField::getRHS() { return F_; }
    Eigen::MatrixXd& PhysicsField::getSolution() { return U_; }

    const Eigen::SparseMatrix<double>& PhysicsField::getStiffnessMatrix() const { return K_; }
    const Eigen::SparseMatrix<double>& PhysicsField::getMassMatrix() const { return M_; }
    const Eigen::MatrixXd& PhysicsField::getRHS() const { return F_; }
    const Eigen::MatrixXd& PhysicsField::getSolution() const { return U_; }

} // namespace Physics
