#ifndef HEAT1D_HPP
#define HEAT1D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>
#include <memory>

namespace Physics {

    class Heat1D : public PhysicsField {
    public:
        // Constructor now takes a material
        explicit Heat1D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

        void setVolumetricHeatSource(const std::vector<double>& source);

    private:
        const Core::Material& material_;
        double k_; // Thermal conductivity, extracted from material
        std::vector<double> volumetric_heat_source_;
    };

} // namespace Physics

#endif // HEAT1D_HPP
