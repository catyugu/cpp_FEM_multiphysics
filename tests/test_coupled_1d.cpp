#include <gtest/gtest.h>
#include <memory>
#include <vector>

// Include all necessary headers from the solver
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/EMag1D.hpp"
#include "physics/Heat1D.hpp"

// Test fixture for our refactored 1D coupled problem
class RefactoredCoupled1DTest : public ::testing::Test {
protected:
    const double length = 1.0;
    const int num_elements = 20;
    const int num_nodes = num_elements + 1;
    const double V0 = 1.0;
    const double VL = 0.0;
    const double T_ambient = 300.0;

    std::unique_ptr<Core::Problem> problem;
    Core::Material copper{"Copper"};

    // SetUp is called before each test
    void SetUp() override {
        // 1. Define material properties
        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);

        // 2. Create the Problem, which handles mesh and DOFManager creation
        // --- FIX IS HERE ---
        // We must first create a named unique_ptr for the mesh to manage its lifetime.
        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));
        // Then, we move ownership of that mesh into the Problem's constructor.
        // This prevents the mesh from being destroyed prematurely.
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        // 3. Add physics fields with the defined material
        problem->addField(std::make_unique<Physics::EMag1D>(copper));
        problem->addField(std::make_unique<Physics::Heat1D>(copper));

        // 4. Finalize the problem setup
        problem->setup();
    }
};

// --- Test Case ---
// This test uses the new 'Problem' class to orchestrate the simulation
// and validates the final temperature distribution against the analytical solution.
TEST_F(RefactoredCoupled1DTest, EndToEndValidationWithProblemClass) {
    // 1. Add Boundary Conditions through the Problem interface
    auto& dof_manager = problem->getDofManager();
    problem->getField("Voltage")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, 0, "Voltage", V0)
    );
    problem->getField("Voltage")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Voltage", VL)
    );
    problem->getField("Temperature")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_ambient)
    );
    problem->getField("Temperature")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Temperature", T_ambient)
    );

    // 2. Solve the entire coupled problem with a single command
    ASSERT_NO_THROW(problem->solve());

    // 3. Validate the final Temperature Solution
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // Calculate analytical solution for max temperature
    double sigma = copper.getProperty("electrical_conductivity");
    double k = copper.getProperty("thermal_conductivity");
    double V_diff = V0 - VL;
    double Q = sigma * (V_diff / length) * (V_diff / length);
    double analytical_max_temp = T_ambient + (Q * length * length) / (8.0 * k);

    // Get the FEM temperature at the center node
    int center_node_id = num_elements / 2;
    int center_eq_idx = dof_manager.getEquationIndex(center_node_id, "Temperature");
    double fem_max_temp = heat_field->getSolution()(center_eq_idx);

    // Assert that the FEM result is close to the analytical one.
    // A tolerance of 1% of the temperature rise is reasonable for this mesh density.
    double tolerance = (analytical_max_temp - T_ambient) * 0.001;
    ASSERT_NEAR(fem_max_temp, analytical_max_temp, tolerance);
}
