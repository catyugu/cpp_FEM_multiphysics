#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat2D.hpp" // Use the new 2D physics

// Test fixture for 2D Heat Conduction
class Heat2DTest : public ::testing::Test {
protected:
    const double size = 1.0; // 1x1 meter plate
    const int num_div = 20;  // 100 divisions along each axis

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("thermal_conductivity", 401.0);

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
    const double T_hot = 500.0;  // K
    const double T_cold = 300.0; // K

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // 1. Apply Boundary Conditions
    // Iterate over all nodes and apply BCs based on position
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        bool is_on_cold_boundary = (std::abs(coords[0] - 0.0) < 1e-9 ||
                                    std::abs(coords[0] - size) < 1e-9 ||
                                    std::abs(coords[1] - 0.0) < 1e-9);

        if (std::abs(coords[1] - size) < 1e-9) { // Top edge is hot
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", T_hot));
        } else if (is_on_cold_boundary) { // Other three edges are cold
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", T_cold));
        }
    }

    // 2. Solve the problem
    ASSERT_NO_THROW(problem->solve());
    problem->exportResults("results_2d_heat.vtk");

    // 3. Validate the solution
    // The analytical solution is complex (a Fourier series).
    // Instead, we will perform plausibility checks based on the Maximum Principle.
    const auto& solution = heat_field->getSolution();

    double max_temp = -1.0;
    double min_temp = 1e10;
    for(int i=0; i<solution.size(); ++i) {
        if(solution(i) > max_temp) max_temp = solution(i);
        if(solution(i) < min_temp) min_temp = solution(i);
    }

    // Check 1: The max/min temperatures should not exceed the boundary values.
    ASSERT_LE(max_temp, T_hot + 1e-9);
    ASSERT_GE(min_temp, T_cold - 1e-9);

    // Check 2: The temperature at the center should be between hot and cold.
    // NOTE: Node ID calculation for the center needs to be robust.
    // Let's find a node near the center instead of guessing its ID.
    Core::Node* center_node = nullptr;
    double min_dist_sq = 1e10;
    for(const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        double dist_sq = std::pow(coords[0] - size/2.0, 2) + std::pow(coords[1] - size/2.0, 2);
        if(dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            center_node = node;
        }
    }
    ASSERT_NE(center_node, nullptr);

    int center_dof = dof_manager.getEquationIndex(center_node->getId(), "Temperature");
    double center_temp = solution(center_dof);
    ASSERT_GT(center_temp, T_cold);
    ASSERT_LT(center_temp, T_hot);

    SimpleLogger::Logger::instance().info("2D Heat test validation passed plausibility checks.");
    SimpleLogger::Logger::instance().info("Center temperature: ", center_temp, " K");
}
