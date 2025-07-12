#include <core/coupling/ElectroThermalCoupling.hpp>
#include "physics/Current1D.hpp"
#include "physics/Heat1D.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "utils/SimpleLogger.hpp"

namespace Core {

    void ElectroThermalCoupling::setup(std::vector<Physics::PhysicsField*>& fields) {
        auto& logger = SimpleLogger::Logger::instance();
        logger.info("Setting up Electro-Thermal coupling...");

        Physics::PhysicsField* emag_field = nullptr;
        Physics::PhysicsField* heat_field = nullptr;

        for (auto& field : fields) {
            if (field->getVariableName() == "Voltage") {
                emag_field = field;
            } else if (field->getVariableName() == "Temperature") {
                heat_field = field;
            }
        }

        if (emag_field && heat_field) {
            if (auto* emag1d = dynamic_cast<Physics::Current1D*>(emag_field)) {
                if (auto* heat1d = dynamic_cast<Physics::Heat1D*>(heat_field)) {
                    emag1d->setCoupledHeatField(heat_field);
                } else {
                    logger.error("Dimensions of the fields must be the same!");
                    return;
                }
            } else if (auto* emag2d = dynamic_cast<Physics::Current2D*>(emag_field)) {
                if (auto* heat2d = dynamic_cast<Physics::Heat2D*>(heat_field)) {
                    emag2d->setCoupledHeatField(heat_field);
                } else {
                    logger.error("Dimensions of the fields must be the same!");
                    return;
                }
            }
        } else {
            logger.error("Could not find both EMag and Heat fields for coupling.");
        }
    }

} // namespace Core