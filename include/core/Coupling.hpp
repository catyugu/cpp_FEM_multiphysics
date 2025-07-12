#ifndef COUPLING_HPP
#define COUPLING_HPP

#include <vector>
#include "physics/PhysicsField.hpp"

namespace Core {

    class Coupling {
    public:
        virtual ~Coupling() = default;
        virtual void setup(std::vector<Physics::PhysicsField*>& fields) = 0;
    };

} // namespace Core

#endif // COUPLING_HPP