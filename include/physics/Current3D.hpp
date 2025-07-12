#ifndef CURRENT3D_HPP
#define CURRENT3D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"

namespace Physics {

    class Current3D : public PhysicsField {
    public:
        explicit Current3D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // CURRENT3D_HPP