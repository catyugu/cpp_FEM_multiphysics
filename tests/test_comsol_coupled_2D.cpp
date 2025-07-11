#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "io/Importer.hpp"
#include "utils/SimpleLogger.hpp"

// Test fixture for the COMSOL benchmark problem
class ComsolBusbarTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    const std::string mesh_filename = "busbar_mesh.mphtxt";

    // --- FIX: Make the Material object a member of the test fixture ---
    // This ensures its lifetime persists for the duration of the test.
    Core::Material copper{"Copper"};

    // This method runs before each test
    void SetUp() override {
        std::unique_ptr<Core::Mesh> mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);
        ASSERT_NE(mesh, nullptr) << "Mesh file '" << mesh_filename << "' could not be found or read.";

        // The material is now a member variable, not a local one.
        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        // Create and set up the Problem
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        // Pass the persistent material object to the constructors
        problem->addField(std::make_unique<Physics::Current2D>(copper));
        problem->addField(std::make_unique<Physics::Heat2D>(copper));

        problem->setup();
    }
};

TEST_F(ComsolBusbarTest, CompareAgainstComsolResult) {
    // --- IMPORTANT ---
    // Replace this value with the maximum temperature you observed in your COMSOL simulation.
    const double COMSOL_MAX_TEMP = 480; // K (Example value, please update!)

    // --- Setup Simulation ---
    const double V_in = 0.1;
    const double T_sink = 293.15;
    const double bar_width = 0.1;

    auto* emag_field = problem->getField("Voltage");
    auto* heat_field = problem->getField("Temperature");
    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();
    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    // --- Apply Boundary Conditions ---
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-4) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", V_in));
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", T_sink));
        } else if (std::abs(coords[0] - bar_width) < 1e-4) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", 0.0));
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", T_sink));
        }
    }

    // --- Solve and Validate ---
    ASSERT_NO_THROW(problem->solveSteadyState());

    // --- Export the results as requested ---
    problem->exportResults("busbar_results.vtk");

    double fem_max_temp = heat_field->getSolution().maxCoeff();
    SimpleLogger::Logger::instance().info("COMSOL Max Temperature: ", COMSOL_MAX_TEMP, " K");
    SimpleLogger::Logger::instance().info("Our FEM Max Temperature: ", fem_max_temp, " K");

    double tolerance = COMSOL_MAX_TEMP * 0.01;
    ASSERT_NEAR(fem_max_temp, COMSOL_MAX_TEMP, tolerance);
}
