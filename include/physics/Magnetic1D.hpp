#ifndef MAGNETIC1D_HPP
#define MAGNETIC1D_HPP

#include "PhysicsField.hpp"
#include "core/Problem.hpp"

namespace Physics {

    class Magnetic1D : public PhysicsField {
    public:
        explicit Magnetic1D();

        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial(const Core::Element* elem) const override {
            return problem_->getMaterial(elem->getMaterialID());
        }
        int getDimension() const override { return 1; }

        const std::string& getFieldVariableName() const override {
            static const std::string var_name = "MagneticVectorPotential";
            return var_name;
        }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

    private:
        Core::Problem* problem_ = nullptr;
    };

} // namespace Physics

#endif // MAGNETIC1D_HPP