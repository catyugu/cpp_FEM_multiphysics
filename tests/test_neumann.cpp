#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"

class Neumann1DTest : public ::testing::Test {
protected:
    const double length = 2.0;
    const int num_elements = 40;
    const int num_nodes = num_elements + 1;
    Core::Material steel{"Steel"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        steel.setProperty("thermal_conductivity", 50.0);
        // --- FIX: Add missing properties ---
        steel.setProperty("density", 7850.0);
        steel.setProperty("specific_heat", 462.0);

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(
           length, num_elements));
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Heat1D>(steel));
        problem->setup();
    }
};

TEST_F(Neumann1DTest, HeatConductionWithFlux) {
    const double T_fixed = 373.15;
    const double flux_out = 1000;

    auto& dof_manager = problem->getDofManager();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_fixed));
    heat_field->addBC(std::make_unique<Core::NeumannBC>(dof_manager, num_nodes - 1, "Temperature", -flux_out));

    ASSERT_NO_THROW(problem->solveSteadyState());

    double k = steel.getProperty("thermal_conductivity");
    const auto& solution_vector = heat_field->getSolution();

    for (int i = 0; i < num_nodes; ++i) {
        double x = problem->getMesh().getNode(i)->getCoords()[0];
        double analytical_T = T_fixed - (flux_out / k) * x;
        double fem_T = solution_vector(dof_manager.getEquationIndex(i, "Temperature"));
        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
}
