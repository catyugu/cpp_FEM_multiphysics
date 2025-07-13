#ifndef EMAG1D_HPP
#define EMAG1D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>

namespace Physics {

    class Current1D : public PhysicsField {
    public:
        explicit Current1D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }
        int getDimension() const override { return 1; }


        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // EMAG1D_HPP