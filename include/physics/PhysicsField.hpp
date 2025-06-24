#ifndef PHYSICSFIELD_HPP
#define PHYSICSFIELD_HPP

#include "core/Mesh.hpp"
#include <vector>

namespace Physics {

    // Abstract base class representing a physical field (e.g., Heat, EM)
    class PhysicsField {
    public:
        virtual ~PhysicsField() = default;

        // Get the name of the physics field
        virtual const char* getName() const = 0;

        // Setup the field for a given mesh
        virtual void setup(Core::Mesh& mesh) = 0;

        // Assemble the global matrices and vectors for the field
        virtual void assemble() = 0;
    };

} // namespace Physics

#endif // PHYSICSFIELD_HPP
