#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Current2D.hpp"
#include "io/Importer.hpp"
#include "utils/SimpleLogger.hpp"

// Test fixture to test the EMag2D field in isolation
class EMag2DSingleFieldTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    const std::string mesh_filename = "busbar_mesh.mphtxt";

    // --- FIX: Make the Material object a member of the test fixture ---
    // This ensures its lifetime persists for the duration of the test.
    Core::Material copper{"Copper"};

    void SetUp() override {
        std::unique_ptr<Core::Mesh> mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);
        ASSERT_NE(mesh, nullptr) << "Mesh file '" << mesh_filename << "' could not be found or read.";

        // The material is now a member variable, not a local one.
        copper.setProperty("electrical_conductivity", 5.96e7);

        problem = std::make_unique<Core::Problem>(std::move(mesh));
        // Pass the persistent material object to the constructor
        problem->addField(std::make_unique<Physics::Current2D>(copper));
        problem->setup();
    }
};

TEST_F(EMag2DSingleFieldTest, SolvesOnImportedMesh) {
    const double V_in = 0.1;
    const double bar_width = 0.1;

    auto* emag_field = problem->getField("Voltage");
    ASSERT_NE(emag_field, nullptr);
    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();

    // Apply Boundary Conditions
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-4) { // Left edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Voltage", V_in));
        } else if (std::abs(coords[0] - bar_width) < 1e-4) { // Right edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(
                dof_manager, node->getId(), "Voltage", 0.0));
        }
    }

    // Solve the problem
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_emag2d_single.vtk");

    // Validate the solution
    const auto& solution = emag_field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        double x = coords[0];
        double analytical_V = V_in * (1.0 - x / bar_width);
        double fem_V = solution(dof_manager.getEquationIndex(node->getId(), "Voltage"));
        ASSERT_NEAR(fem_V, analytical_V, 1e-9);
    }
    SimpleLogger::Logger::instance().info("EMag2D single field test passed successfully.");
}
