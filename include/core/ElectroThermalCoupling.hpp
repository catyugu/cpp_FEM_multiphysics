#ifndef ELECTROTHERMALCOUPLING_HPP
#define ELECTROTHERMALCOUPLING_HPP

#include "Coupling.hpp"

namespace Core {

    class ElectroThermalCoupling : public Coupling {
    public:
        void setup(std::vector<Physics::PhysicsField*>& fields) override;
    };

} // namespace Core

#endif // ELECTROTHERMALCOUPLING_HPP