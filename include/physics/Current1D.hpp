#ifndef EMAG1D_HPP
#define EMAG1D_HPP

#include "PhysicsField.hpp"
#include "core/Problem.hpp" // Include Problem to access materials
#include <vector>

namespace Physics {

    class Current1D : public PhysicsField {
    public:
        explicit Current1D();

        const char* getName() const override;
        const char* getVariableName() const override;
        int getDimension() const override { return 1; }
        const Core::Material& getMaterial(const Core::Element* elem) const override {
            return problem_->getMaterial(elem->getMaterialID());
        }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble(const PhysicsField *coupled_field) override;
    };

} // namespace Physics

#endif // EMAG1D_HPP