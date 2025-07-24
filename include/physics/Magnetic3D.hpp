#ifndef MAGNETIC3D_HPP
#define MAGNETIC3D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"

namespace Physics {

    class Magnetic3D : public PhysicsField {
    public:
        explicit Magnetic3D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }
        int getDimension() const override { return 3; }
        int getNumComponents() const override { return 3; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // MAGNETIC3D_HPP