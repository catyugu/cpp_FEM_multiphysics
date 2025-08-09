#ifndef HEAT2D_HPP
#define HEAT2D_HPP

#include "physics/PhysicsField.hpp"
#include "core/Problem.hpp"
#include <vector>

namespace Physics {

    class Heat2D : public PhysicsField {
    public:
        explicit Heat2D();

        const char* getName() const override;
        const char* getVariableName() const override;
        int getDimension() const override { return 2; }

        const std::string& getFieldVariableName() const override {
            static const std::string var_name = "Temperature";
            return var_name;
        }

        void setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) override;
        void assemble() override;

        const Core::Material& getMaterial(const Core::Element* elem) const override {
            return problem_->getMaterial(elem->getMaterialID());
        }
    };

} // namespace Physics

#endif // HEAT2D_HPP