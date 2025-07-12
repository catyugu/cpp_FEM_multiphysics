#ifndef HEAT3D_HPP
#define HEAT3D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"

namespace Physics {

    class Heat3D : public PhysicsField {
    public:
        explicit Heat3D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
        double k_; // Isotropic thermal conductivity
    };

} // namespace Physics

#endif // HEAT3D_HPP
