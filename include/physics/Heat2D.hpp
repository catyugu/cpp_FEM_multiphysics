#ifndef HEAT2D_HPP
#define HEAT2D_HPP

#include "physics/PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>

namespace Physics {

    class Heat2D : public PhysicsField {
    public:
        explicit Heat2D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }
        int getDimension() const override { return 2; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble(const PhysicsField *coupled_field) override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // HEAT2D_HPP