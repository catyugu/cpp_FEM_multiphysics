#ifndef EMAG2D_HPP
#define EMAG2D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>

namespace Physics {

    class Current2D : public PhysicsField {
    public:
        explicit Current2D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }
        int getDimension() const override { return 2; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // EMAG2D_HPP