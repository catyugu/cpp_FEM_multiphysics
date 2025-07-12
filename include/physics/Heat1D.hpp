#ifndef HEAT1D_HPP
#define HEAT1D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>

namespace Physics {

    class Heat1D : public PhysicsField {
    public:
        explicit Heat1D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // HEAT1D_HPP