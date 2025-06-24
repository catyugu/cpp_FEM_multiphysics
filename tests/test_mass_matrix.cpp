#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Heat2D.hpp"

// Test fixture for mass matrix validation
class MassMatrixTest : public ::testing::Test {
protected:
    const double width = 0.5;
    const double height = 0.2;
    const int nx = 2;
    const int ny = 1;

    Core::Material test_mat{"TestMaterial"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        // Use properties that are easy to verify
        test_mat.setProperty("thermal_conductivity", 1.0);
        test_mat.setProperty("density", 100.0);
        test_mat.setProperty("specific_heat", 2.0);

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(width, height, nx, ny)
        );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Heat2D>(test_mat));
        problem->setup();
    }
};

TEST_F(MassMatrixTest, MassMatrixSummation) {
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    // Assemble the K and M matrices
    heat_field->assemble();

    const auto& M = heat_field->getMassMatrix();

    // The sum of all entries in the consistent mass matrix should equal rho * cp * total_volume.
    // For a 2D problem (assuming unit thickness), this is rho * cp * total_area.
    double total_mass = M.sum();

    double rho = test_mat.getProperty("density");
    double cp = test_mat.getProperty("specific_heat");
    double total_area = width * height;

    double expected_total_mass = rho * cp * total_area;

    SimpleLogger::Logger::instance().info("Total mass from matrix sum: ", total_mass);
    SimpleLogger::Logger::instance().info("Expected total mass (rho*cp*A): ", expected_total_mass);

    // Verify that the computed total mass matches the expected value.
    ASSERT_NEAR(total_mass, expected_total_mass, 1e-9);
}
