#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat2D.hpp"

// Test fixture for 2D Heat Conduction
class Heat2DTest : public ::testing::Test {
protected:
    const double size = 1.0;
    const int num_div = 10;

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("thermal_conductivity", 401.0);
        // --- FIX: Add missing properties required by the updated Heat2D class ---
        copper.setProperty("density", 8960.0);       // Value doesn't matter for steady-state
        copper.setProperty("specific_heat", 385.0);  // Value doesn't matter for steady-state

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(size, size, num_div, num_div)
        );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Heat2D>(copper));
        problem->setup();
    }
};

// Test a square plate with fixed temperatures on its boundaries
TEST_F(Heat2DTest, SquarePlateFixedBoundaries) {
    const double T_hot = 500.0;
    const double T_cold = 300.0;

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // Apply Boundary Conditions
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        bool is_on_cold_boundary = (std::abs(coords[0] - 0.0) < 1e-9 ||
                                    std::abs(coords[0] - size) < 1e-9 ||
                                    std::abs(coords[1] - 0.0) < 1e-9);

        if (std::abs(coords[1] - size) < 1e-9) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", T_hot));
        } else if (is_on_cold_boundary) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", T_cold));
        }
    }

    // Solve and validate
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_2d_heat.vtk");

    const auto& solution = heat_field->getSolution();
    double max_temp = -1.0, min_temp = 1e10;
    for(int i=0; i<solution.size(); ++i) {
        if(solution(i) > max_temp) max_temp = solution(i);
        if(solution(i) < min_temp) min_temp = solution(i);
    }

    ASSERT_LE(max_temp, T_hot + 1e-9);
    ASSERT_GE(min_temp, T_cold - 1e-9);
}
