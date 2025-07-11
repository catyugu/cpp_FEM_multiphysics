#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat3D.hpp" // Use the new 3D physics

// Test fixture for 3D Heat Conduction
class Heat3DTest : public ::testing::Test {
protected:
    const double size = 1.0; // 1x1x1 meter cube
    const int num_div = 5;   // 5 divisions along each axis

    Core::Material steel{"Steel"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        steel.setProperty("thermal_conductivity", 50.2);
        steel.setProperty("density", 7850.0);
        steel.setProperty("specific_heat", 462.0);
        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_3d_mesh(size, size, size, num_div, num_div, num_div));
        problem = std::make_unique<Core::Problem>(
           std::move(mesh)
        );

        problem->addField(std::make_unique<Physics::Heat3D>(steel));
        problem->setup();
    }
};

// Test a cube with a fixed temperature gradient across the x-axis
TEST_F(Heat3DTest, CubeFixedBoundaries) {
    const double T_hot = 400.0;  // K
    const double T_cold = 300.0; // K

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // 1. Apply Boundary Conditions
    // T_hot on the x=0 face, T_cold on the x=1 face
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) { // Left face
            heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_hot)));
        } else if (std::abs(coords[0] - size) < 1e-9) { // Right face
             heat_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_cold)));
        }
    }

    // 2. Solve the problem
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_3d_heat.vtk");

    // 3. Validate the solution against the linear analytical solution
    // T(x) = T_hot - (T_hot - T_cold) * (x / size)
    const auto& solution = heat_field->getSolution();
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        double x = coords[0];

        // Skip validation for boundary nodes as their values are fixed
        if (std::abs(x - 0.0) < 1e-9 || std::abs(x - size) < 1e-9) {
            continue;
        }

        double analytical_T = T_hot - (T_hot - T_cold) * (x / size);
        int dof_idx = dof_manager.getEquationIndex(node->getId(), "Temperature");
        double fem_T = solution(dof_idx);

        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
    SimpleLogger::Logger::instance().info("3D Heat test passed. FEM solution matches analytical solution.");
}
