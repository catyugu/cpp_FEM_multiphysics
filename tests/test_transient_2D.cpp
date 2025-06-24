#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Heat2D.hpp"
#include "core/BoundaryCondition.hpp"
#include "utils/SimpleLogger.hpp"

// Test fixture for 2D transient heat conduction
class Transient2DTest : public ::testing::Test {
protected:
    const double size = 1.0;
    const int num_div = 10;

    Core::Material aluminum{"Aluminum"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        aluminum.setProperty("thermal_conductivity", 237.0);
        aluminum.setProperty("density", 2700.0);
        aluminum.setProperty("specific_heat", 900.0);
        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(size, size, num_div, num_div)
            );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Heat2D>(aluminum));
        problem->setup();
    }
};

// This test simulates the cooling of a hot plate and checks that the center cools down.
TEST_F(Transient2DTest, PlateCooling) {
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    const double T_initial = 600.0; // K
    const double T_boundary = 300.0; // K

    // 1. Set simulation parameters for a short duration
    problem->setTimeSteppingParameters(0.1, 5.0); // 50 steps of 0.1s each

    // 2. Set initial and boundary conditions
    heat_field->setInitialConditions(T_initial);

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();

    // Apply fixed temperature to all boundary nodes
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        bool is_boundary = (std::abs(coords[0] - 0.0) < 1e-9 ||
                            std::abs(coords[0] - size) < 1e-9 ||
                            std::abs(coords[1] - 0.0) < 1e-9 ||
                            std::abs(coords[1] - size) < 1e-9);
        if (is_boundary) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", T_boundary));
        }
    }

    // 3. Solve the transient problem
    problem->solveTransient();
    problem->exportResults("results_2d_transient.vtk");

    // 4. Validate the result
    const auto& solution = heat_field->getSolution();

    // Find the center node
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
    double final_center_temp = solution(center_dof);

    SimpleLogger::Logger::instance().info("Initial center temperature: ", T_initial, " K");
    SimpleLogger::Logger::instance().info("Final center temperature:   ", final_center_temp, " K");
    SimpleLogger::Logger::instance().info("Boundary temperature:       ", T_boundary, " K");

    // The final temperature at the center must be lower than the initial temperature
    // but higher than the boundary temperature, as it has not reached steady state.
    ASSERT_LT(final_center_temp, T_initial);
    ASSERT_GT(final_center_temp, T_boundary);
}
