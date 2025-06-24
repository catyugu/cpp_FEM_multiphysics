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

        // Sets the volumetric heat source for each element.
        void setVolumetricHeatSource(const std::vector<double>& source);

    private:
        const Core::Material& material_;
        double k_; // Isotropic thermal conductivity
        std::vector<double> volumetric_heat_source_;
    };

} // namespace Physics

#endif // HEAT2D_HPP
