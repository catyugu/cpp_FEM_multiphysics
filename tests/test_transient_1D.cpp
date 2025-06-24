#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Heat1D.hpp"
#include "core/BoundaryCondition.hpp"
#include "utils/SimpleLogger.hpp" // Include logger for output

class Transient1DTest : public ::testing::Test {
protected:
    const double length = 1.0;
    const int num_elements = 20;

    Core::Material steel{"Steel"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        steel.setProperty("thermal_conductivity", 50.2);
        steel.setProperty("density", 7850.0);
        steel.setProperty("specific_heat", 462.0);

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));
        problem =  std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Heat1D>(steel));
        problem->setup();
    }
};

// This test verifies that the transient solution approaches the correct
// steady-state solution over a long period.
TEST_F(Transient1DTest, ConvergesToSteadyState) {
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // 1. Set simulation parameters
    problem->setTimeSteppingParameters(50.0, 100000.0); // Large time steps, long total time

    // 2. Set initial and boundary conditions
    heat_field->setInitialConditions(0.0); // Start at 0 Kelvin

    auto& dof_manager = problem->getDofManager();
    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", 100.0));
    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_elements, "Temperature", 0.0));

    // 3. Solve the transient problem using the generic solver
    problem->solveTransient();

    // 4. Validate the final state against the analytical steady-state solution
    // T(x) = 100 * (1 - x/L)
    const auto& solution = heat_field->getSolution();
    const auto& mesh = problem->getMesh();
    for(int i = 0; i <= num_elements; ++i) {
        double x = mesh.getNode(i)->getCoords()[0];
        double analytical_T = 100.0 * (1.0 - x / length);
        double fem_T = solution(dof_manager.getEquationIndex(i, "Temperature"));

        ASSERT_NEAR(fem_T, analytical_T, 1e-2); // Check with a reasonable tolerance
    }
     SimpleLogger::Logger::instance().info("Transient test passed: Converged to correct steady state.");
}
