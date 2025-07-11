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

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

        std::vector<double> calculateJouleHeat() const;

        // Links this EMag field to a heat field for temperature-dependent calculations.
        void setCoupledHeatField(const PhysicsField* heat_field);

    private:
        const Core::Material& material_;
        const PhysicsField* heat_field_ = nullptr; // FIX: Add the missing member declaration
    };

} // namespace Physics

#endif // EMAG1D_HPP
