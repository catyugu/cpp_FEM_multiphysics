#ifndef EMAG2D_HPP
#define EMAG2D_HPP

#include "PhysicsField.hpp"
#include "core/Problem.hpp"
#include <vector>

namespace Physics {

    class Current2D : public PhysicsField {
    public:
        explicit Current2D();

        const char* getName() const override;
        const char* getVariableName() const override;
        int getDimension() const override { return 2; }
        const Core::Material& getMaterial(const Core::Element* elem) const override {
            return problem_->getMaterial(elem->getMaterialID());
        }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble(const PhysicsField *coupled_field) override;

    private:
    };

} // namespace Physics

#endif // EMAG2D_HPP