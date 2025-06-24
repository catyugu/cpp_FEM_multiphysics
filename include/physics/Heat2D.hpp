#ifndef HEAT2D_HPP
#define HEAT2D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>

namespace Physics {

    class Heat2D : public PhysicsField {
    public:
        explicit Heat2D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

        void setVolumetricHeatSource(const std::vector<double>& source);

    private:
        const Core::Material& material_;
        double k_;     // Thermal conductivity
        double rho_;   // Density
        double cp_;    // Specific heat capacity
        std::vector<double> volumetric_heat_source_;
    };

} // namespace Physics

#endif // HEAT2D_HPP
