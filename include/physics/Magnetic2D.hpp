#ifndef MAGNETIC2D_HPP
#define MAGNETIC2D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"

namespace Physics {

    class Magnetic2D : public PhysicsField {
    public:
        explicit Magnetic2D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // MAGNETIC2D_HPP