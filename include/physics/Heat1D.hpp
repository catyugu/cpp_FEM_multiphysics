#ifndef HEAT1D_HPP
#define HEAT1D_HPP

#include "PhysicsField.hpp"

namespace Physics {

    class Heat1D : public PhysicsField {
    public:
        Heat1D(double k);

        const char* getName() const override;
        void setup(Core::Mesh& mesh) override;
        void assemble() override;

        // Method to set the heat source from another field (coupling)
        void setHeatSource(const std::vector<double>& source);


    private:
        Core::Mesh* mesh_;
        double k_; // Thermal conductivity
        std::vector<double> heat_source_;
    };

} // namespace Physics

#endif // HEAT1D_HPP
