#ifndef ELECTROTHERMALCOUPLING_HPP
#define ELECTROTHERMALCOUPLING_HPP

#include "Coupling.hpp"

namespace Physics {
    class Current2D;
    class Heat2D;
}
namespace Core {
    class ElectroThermalCoupling : public Coupling {
    public:
        void setup(std::vector<Physics::PhysicsField*>& fields) override;
        void execute() override;

    private:
        Physics::PhysicsField* emag_field_ = nullptr;
        Physics::PhysicsField* heat_field_ = nullptr;
    };
}
#endif // ELECTROTHERMALCOUPLING_HPP