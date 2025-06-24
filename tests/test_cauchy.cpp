#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat1D.hpp"

// Test fixture for Cauchy (mixed) boundary condition validation
class Cauchy1DTest : public ::testing::Test {
protected:
    const double length = 2.0;    // m
    const int num_elements = 40;
    const int num_nodes = num_elements + 1;

    // Define a simple material
    Core::Material aluminum{"Aluminum"};

    // Test parameters
    const double T_fixed = 400.0;     // K (fixed temperature at x=0)
    const double T_ambient = 293.15;  // K (ambient temperature for convection)
    const double h_conv = 15.0;       // W/(m^2*K) (convection coefficient)

    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        aluminum.setProperty("thermal_conductivity", 237.0); // W/(m*K)

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));

        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Heat1D>(aluminum));
        problem->setup();
    }
};

// Test a 1D heat conduction problem with one Dirichlet and one Cauchy BC
TEST_F(Cauchy1DTest, HeatConductionWithConvection) {
    auto& dof_manager = problem->getDofManager();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // 1. Apply Boundary Conditions
    // Fixed temperature T=T_fixed at x=0
    heat_field->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_fixed)
    );
    // Convection (Cauchy) condition at x=L
    heat_field->addBC(
        std::make_unique<Core::CauchyBC>(dof_manager, num_nodes - 1, "Temperature", h_conv, T_ambient)
    );

    // 2. Solve the problem
    ASSERT_NO_THROW(problem->solve());

    // 3. Validate the solution against the analytical result
    // Analytical solution: T(x) = C1*x + C2
    // C2 = T_fixed
    // C1 = h*(T_ambient - T_fixed) / (h*L + k)
    double k = aluminum.getProperty("thermal_conductivity");
    double C1 = (h_conv * (T_ambient - T_fixed)) / (h_conv * length + k);
    double C2 = T_fixed;

    const auto& solution_vector = heat_field->getSolution();

    SimpleLogger::Logger::instance().info("\n--- Cauchy BC Test Validation ---");
    for (int i = 0; i < num_nodes; ++i) {
        double x = problem->getMesh().getNode(i)->getCoords()[0];
        double analytical_T = C1 * x + C2;

        int dof_idx = dof_manager.getEquationIndex(i, "Temperature");
        double fem_T = solution_vector(dof_idx);

        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
    SimpleLogger::Logger::instance().info("Cauchy BC test passed. FEM solution matches analytical solution.");
}
