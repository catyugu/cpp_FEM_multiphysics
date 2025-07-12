#ifndef COUPLING_HPP
#define COUPLING_HPP

#include <vector>
#include "physics/PhysicsField.hpp"

namespace Core {

    class Coupling {
    public:
        virtual ~Coupling() = default;
        virtual void setup(std::vector<Physics::PhysicsField*>& fields) = 0;
        virtual void execute() = 0; // New method
    };

} // namespace Core

#endif // COUPLING_HPP