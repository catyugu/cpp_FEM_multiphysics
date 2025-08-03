#ifndef HEAT1D_HPP
#define HEAT1D_HPP

#include "PhysicsField.hpp"
#include "core/Problem.hpp"
#include <vector>

namespace Physics {

    class Heat1D : public PhysicsField {
    public:
        explicit Heat1D();
        
        const char* getName() const override;
        const char* getVariableName() const override;
        const Core::Material& getMaterial(const Core::Element* elem) const override {
            return problem_->getMaterial(elem->getMaterialID());
        }
        int getDimension() const override { return 1; }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble(const PhysicsField *coupled_field) override;
    };

} // namespace Physics

#endif // HEAT1D_HPP