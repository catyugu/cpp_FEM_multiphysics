#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"

// Test fixture for 2D Coupled Electro-Thermal simulation
class Coupled2DTest : public ::testing::Test {
protected:
    const double width = 0.02;
    const double height = 0.01;
    const int nx = 20;
    const int ny = 10;

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);
        // --- FIX: Add missing properties required by the updated Heat2D class ---
        copper.setProperty("density", 8960.0);       // Value doesn't affect steady-state result
        copper.setProperty("specific_heat", 385.0);  // Value doesn't affect steady-state result

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(width, height, nx, ny)
        );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Current2D>(copper));
        problem->addField(std::make_unique<Physics::Heat2D>(copper));

        problem->setup();
    }
};

// Test a rectangular plate with Joule heating and fixed temperature boundaries
TEST_F(Coupled2DTest, RectangularPlateJouleHeating) {
    const double V_high = 1.0, V_low = 0.0;
    const double T_boundary = 300.0;

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    auto* emag_field = problem->getField("Voltage");
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    // Apply BCs
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        bool is_boundary = (std::abs(coords[0] - 0.0) < 1e-9 ||
                            std::abs(coords[0] - width) < 1e-9 ||
                            std::abs(coords[1] - 0.0) < 1e-9 ||
                            std::abs(coords[1] - height) < 1e-9);

        if (std::abs(coords[0] - 0.0) < 1e-9) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_high)));
        } else if (std::abs(coords[0] - width) < 1e-9) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_low)));
        }

        if (is_boundary) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_boundary)));
        }
    }

    // Solve and validate
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_2d_coupled.vtk");

    const auto& temp_solution = heat_field->getSolution();
    double max_temp = temp_solution.maxCoeff();

    ASSERT_GT(max_temp, T_boundary);
}
