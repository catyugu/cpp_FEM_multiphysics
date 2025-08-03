#ifndef MAGNETIC2D_HPP
#define MAGNETIC2D_HPP

#include "PhysicsField.hpp"
#include "core/Problem.hpp"

namespace Physics {

    class Magnetic2D : public PhysicsField {
    public:
        explicit Magnetic2D();
        
        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial(const Core::Element* elem) const override;
        int getDimension() const override { return 2; }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble(const PhysicsField *coupled_field) override;
    };

} // namespace Physics

#endif // MAGNETIC2D_HPP