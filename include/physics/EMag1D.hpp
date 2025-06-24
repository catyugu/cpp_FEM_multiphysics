#ifndef EMAG1D_HPP
#define EMAG1D_HPP

#include "PhysicsField.hpp"
#include "core/Material.hpp"
#include <vector>
#include <memory>

namespace Physics {

    class EMag1D : public PhysicsField {
    public:
        // Constructor now takes a material
        explicit EMag1D(const Core::Material& material);

        const char* getName() const override;
        const char* getVariableName() const override;

        void setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

        std::vector<double> calculateJouleHeat() const;

    private:
        const Core::Material& material_;
        double sigma_; // Electrical conductivity, extracted from material
    };

} // namespace Physics

#endif // EMAG1D_HPP
