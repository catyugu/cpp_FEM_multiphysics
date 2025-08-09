#ifndef MAGNETIC3D_HPP
#define MAGNETIC3D_HPP

#include "PhysicsField.hpp"
#include "core/Problem.hpp"

namespace Physics {

    class Magnetic3D : public PhysicsField {
    public:
        explicit Magnetic3D();

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial(const Core::Element* elem) const override {
            return problem_->getMaterial(elem->getMaterialID());
        }
        int getDimension() const override { return 3; }
        int getNumComponents() const override { return 3; }

        const std::string& getFieldVariableName() const override {
            static const std::string var_name = "MagneticVectorPotential";
            return var_name;
        }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        // No material_ member anymore, fetched via Problem
    };

} // namespace Physics

#endif // MAGNETIC3D_HPP