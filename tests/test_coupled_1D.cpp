#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <core/LinearSolver.hpp>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/EMag1D.hpp"
#include "physics/Heat1D.hpp"

// Test fixture for 1D coupled problem
class Coupled1DTest : public ::testing::Test {
protected:
    const double length = 1.0;
    const int num_elements = 20;
    const int num_nodes = num_elements + 1;

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);
        // --- FIX: Add missing properties ---
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(
            length, num_elements));
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::EMag1D>(copper));
        problem->addField(std::make_unique<Physics::Heat1D>(copper));

        problem->setup();
    }
};

TEST_F(Coupled1DTest, EndToEndValidation) {
    const double V0 = 1.0, VL = 0.0;
    const double T_ambient = 300.0;

    auto& dof_manager = problem->getDofManager();
    problem->getField("Voltage")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Voltage", V0));
    problem->getField("Voltage")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Voltage", VL));
    problem->getField("Temperature")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_ambient));
    problem->getField("Temperature")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Temperature", T_ambient));

    // For a 1D coupled problem, we need a custom iterative solve for now
    // as Problem::solveSteadyState is geared towards 2D.
    // This logic can be moved into the Problem class later.
    auto* emag_field = problem->getField("Voltage");
    auto* heat_field = problem->getField("Temperature");

    emag_field->assemble();
    emag_field->applyBCs();
    Core::LinearSolver::solve(emag_field->getStiffnessMatrix(), emag_field->getRHSVector(), emag_field->getSolution());

    auto joule_heat = dynamic_cast<Physics::EMag1D*>(emag_field)->calculateJouleHeat();
    dynamic_cast<Physics::Heat1D*>(heat_field)->setVolumetricHeatSource(joule_heat);

    heat_field->assemble();
    heat_field->applyBCs();
    Core::LinearSolver::solve(heat_field->getStiffnessMatrix(), heat_field->getRHSVector(), heat_field->getSolution());

    double k = copper.getProperty("thermal_conductivity");
    double analytical_max_temp = T_ambient + (copper.getProperty("electrical_conductivity") / (8.0 * k)) * std::pow(V0 - VL, 2);
    double fem_max_temp = heat_field->getSolution().maxCoeff();

    ASSERT_NEAR(fem_max_temp, analytical_max_temp, analytical_max_temp * 0.01);
}
