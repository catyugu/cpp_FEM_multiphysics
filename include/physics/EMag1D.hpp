#ifndef EMAG1D_HPP
#define EMAG1D_HPP

#include "PhysicsField.hpp"

namespace Physics {

    class EMag1D : public PhysicsField {
    public:
        EMag1D(double sigma);

        const char* getName() const override;
        void setup(Core::Mesh& mesh) override;
        void assemble() override;

    private:
        Core::Mesh* mesh_;
        double sigma_; // Electrical conductivity
        // Here you would have your global stiffness matrix and RHS vector
        // For now, we'll just log the assembly process.
    };

} // namespace Physics

#endif // EMAG1D_HPP
