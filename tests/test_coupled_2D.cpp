#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/EMag2D.hpp"
#include "physics/Heat2D.hpp"

// Test fixture for 2D Coupled Electro-Thermal simulation
class Coupled2DTest : public ::testing::Test {
protected:
    const double width = 0.02; // 2 cm
    const double height = 0.01; // 1 cm
    const int nx = 20;
    const int ny = 10;

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(width, height, nx, ny)
        );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        // Add both 2D physics fields
        problem->addField(std::make_unique<Physics::EMag2D>(copper));
        problem->addField(std::make_unique<Physics::Heat2D>(copper));

        problem->setup();
    }
};

// Test a rectangular plate with Joule heating and fixed temperature boundaries
TEST_F(Coupled2DTest, RectangularPlateJouleHeating) {
    const double V_high = 1.0;  // Volts
    const double V_low = 0.0;   // Volts
    const double T_boundary = 300.0; // K

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    auto* emag_field = problem->getField("Voltage");
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    // 1. Apply Boundary Conditions
    // Voltage difference on left/right edges, fixed temperature on all four edges.
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        bool is_boundary = (std::abs(coords[0] - 0.0) < 1e-9 ||
                            std::abs(coords[0] - width) < 1e-9 ||
                            std::abs(coords[1] - 0.0) < 1e-9 ||
                            std::abs(coords[1] - height) < 1e-9);

        // Voltage BCs
        if (std::abs(coords[0] - 0.0) < 1e-9) { // Left edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Voltage", V_high));
        } else if (std::abs(coords[0] - width) < 1e-9) { // Right edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Voltage", V_low));
        }

        // Temperature BCs
        if (is_boundary) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", T_boundary));
        }
    }

    // 2. Solve the problem
    ASSERT_NO_THROW(problem->solve());
    problem->exportResults("results_2d_coupled.vtk");

    // 3. Validate the solution using plausibility checks
    const auto& temp_solution = heat_field->getSolution();

    double max_temp = 0.0;
    for(int i=0; i<temp_solution.size(); ++i) {
        if(temp_solution(i) > max_temp) {
            max_temp = temp_solution(i);
        }
    }

    // Plausibility Check 1: The maximum temperature must be greater than the boundary temp
    // because heat is generated internally.
    ASSERT_GT(max_temp, T_boundary);

    // Plausibility Check 2: The maximum temperature should occur roughly in the center.
    Core::Node* center_node = nullptr;
    double min_dist_sq = 1e10;
    for(const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        double dist_sq = std::pow(coords[0] - width/2.0, 2) + std::pow(coords[1] - height/2.0, 2);
        if(dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            center_node = node;
        }
    }
    ASSERT_NE(center_node, nullptr);
    int center_dof = dof_manager.getEquationIndex(center_node->getId(), "Temperature");
    double center_temp = temp_solution(center_dof);

    // The center temperature should be very close to the maximum temperature in the domain.
    ASSERT_NEAR(center_temp, max_temp, max_temp * 0.05); // Allow 5% tolerance

    SimpleLogger::Logger::instance().info("2D Coupled test validation passed plausibility checks.");
    SimpleLogger::Logger::instance().info("Max temperature: ", max_temp, " K");
}
