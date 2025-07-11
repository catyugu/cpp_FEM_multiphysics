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

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

        // Calculates the volumetric Joule heat for each element.
        std::vector<double> calculateJouleHeat() const;

        /**
         * @brief Links this EMag field to a heat field to enable temperature-dependent calculations.
         * @param heat_field A pointer to the problem's heat field.
         */
        void setCoupledHeatField(const PhysicsField* heat_field);

    private:
        const Core::Material& material_;
        const PhysicsField* heat_field_ = nullptr; // Pointer to the heat field for temperature data
    };

} // namespace Physics

#endif // EMAG2D_HPP
