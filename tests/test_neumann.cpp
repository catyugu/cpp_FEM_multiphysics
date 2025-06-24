#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat1D.hpp"

// Test fixture for Neumann boundary condition validation
class Neumann1DTest : public ::testing::Test {
protected:
    const double length = 2.0;
    const int num_elements = 40;
    const int num_nodes = num_elements + 1;

    // Define a simple material
    Core::Material steel{"Steel"};

    // Test parameters
    const double T_fixed = 373.15; // 100 C
    const double flux_out = 1000;  // W/m^2 (heat leaving the rod)

    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        steel.setProperty("thermal_conductivity", 50.0); // W/(m*K)

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));

        problem = std::make_unique<Core::Problem>(std::move(mesh));


        // Add only the heat transfer physics field
        problem->addField(std::make_unique<Physics::Heat1D>(steel));

        problem->setup();
    }
};

// Test a 1D heat conduction problem with one Dirichlet and one Neumann BC
TEST_F(Neumann1DTest, HeatConductionWithFlux) {
    auto& dof_manager = problem->getDofManager();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // 1. Apply Boundary Conditions
    // Fixed temperature T=T_fixed at x=0
    heat_field->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_fixed)
    );
    // Outward flux q=-1000 at x=L. In our formulation, this corresponds to a
    // negative value in the force vector.
    heat_field->addBC(
        std::make_unique<Core::NeumannBC>(dof_manager, num_nodes - 1, "Temperature", -flux_out)
    );

    // 2. Solve the problem
    ASSERT_NO_THROW(problem->solve());

    // 3. Validate the solution against the analytical result
    // Analytical solution for this case: T(x) = T_fixed - (flux_out / k) * x
    double k = steel.getProperty("thermal_conductivity");
    const auto& solution_vector = heat_field->getSolution();

    SimpleLogger::Logger::instance().info("\n--- Neumann BC Test Validation ---");
    for (int i = 0; i < num_nodes; ++i) {
        double x = problem->getMesh().getNode(i)->getCoords()[0];
        double analytical_T = T_fixed - (flux_out / k) * x;

        int dof_idx = dof_manager.getEquationIndex(i, "Temperature");
        double fem_T = solution_vector(dof_idx);

        // Use a small absolute tolerance for floating point comparisons
        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
    SimpleLogger::Logger::instance().info("Neumann BC test passed. FEM solution matches analytical solution.");
}
