#ifndef COUPLINGMANAGER_HPP
#define COUPLINGMANAGER_HPP

#include <vector>
#include <memory>
#include "physics/PhysicsField.hpp"
#include "Coupling.hpp"

namespace Core {

    class CouplingManager {
    public:
        void registerField(Physics::PhysicsField& field);
        void addCoupling(std::unique_ptr<Coupling> coupling);
        void setupCouplings();

    private:
        std::vector<Physics::PhysicsField*> fields_;
        std::vector<std::unique_ptr<Coupling>> couplings_;
    };

} // namespace Core

#endif // COUPLINGMANAGER_HPP