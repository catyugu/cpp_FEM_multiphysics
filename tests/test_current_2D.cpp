#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Current2D.hpp" // Use the new 2D physics class

// Test fixture for 2D Electrostatics
class Current2DTest : public ::testing::Test {
protected:
    const double width = 2.0;
    const double height = 1.0;
    const int nx = 20;
    const int ny = 10;

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("electrical_conductivity", 5.96e7);

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(width, height, nx, ny)
            );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Current2D>(copper));
        problem->setup();
    }
};

// Test a rectangular plate with a fixed voltage difference
TEST_F(Current2DTest, RectangularPlateVoltageDrop) {
    const double V_high = 5.0; // Volts
    const double V_low = 1.0;  // Volts

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    auto* emag_field = problem->getField("Voltage");
    ASSERT_NE(emag_field, nullptr);

    // 1. Apply Boundary Conditions
    // Set V_high on the left edge (x=0) and V_low on the right edge (x=width)
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) { // Left edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_high)));
        } else if (std::abs(coords[0] - width) < 1e-9) { // Right edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_low)));
        }
    }

    // 2. Solve the problem
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_2d_emag.vtk");

    // 3. Validate the solution against the analytical result
    // The analytical solution is a linear potential drop: V(x) = V_high - (V_high - V_low) * (x / width)
    const auto& solution = emag_field->getSolution();

    SimpleLogger::Logger::instance().info("\n--- 2D EMag Test Validation ---");
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        double x = coords[0];
        double analytical_V = V_high - (V_high - V_low) * (x / width);

        int dof_idx = dof_manager.getEquationIndex(node->getId(), "Voltage");
        double fem_V = solution(dof_idx);

        ASSERT_NEAR(fem_V, analytical_V, 1e-9);
    }
     SimpleLogger::Logger::instance().info("2D EMag test passed. FEM solution matches analytical solution.");
}
