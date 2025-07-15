#include <core/coupling/CouplingManager.hpp>
#include "utils/SimpleLogger.hpp"

namespace Core {

    void CouplingManager::registerField(Physics::PhysicsField& field) {
        fields_.push_back(&field);
    }

    void CouplingManager::addCoupling(std::unique_ptr<Coupling> coupling) {
        couplings_.push_back(std::move(coupling));
    }

    void CouplingManager::setupCouplings() {
        auto& logger = Utils::Logger::instance();
        logger.info("Setting up couplings...");
        for (const auto& coupling : couplings_) {
            coupling->setup(fields_);
        }
        logger.info("Couplings setup complete.");
    }

    void CouplingManager::executeCouplings() {
        for (const auto& coupling : couplings_) {
            coupling->execute();
        }
    }

} // namespace Core