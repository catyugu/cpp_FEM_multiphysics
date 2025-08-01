#ifndef MAGNETIC1D_HPP
#define MAGNETIC1D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"

namespace Physics {

    class Magnetic1D : public PhysicsField {
    public:
        explicit Magnetic1D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial() const override { return material_; }
        int getDimension() const override { return 1; }

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble(const PhysicsField *coupled_field) override;

    private:
        const Core::Material& material_;
    };

} // namespace Physics

#endif // MAGNETIC1D_HPP