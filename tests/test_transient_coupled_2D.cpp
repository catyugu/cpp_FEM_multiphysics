#include <gtest/gtest.h>
#include <memory>
#include <solver/SolverFactory.hpp>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "core/BoundaryCondition.hpp"
#include "utils/SimpleLogger.hpp"

// Test fixture for 2D transient coupled simulation
class TransientCoupled2DTest : public ::testing::Test {
protected:
    const double width = 0.02, height = 0.01;
    const int nx = 20, ny = 10;

    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        // Using temperature-dependent model for conductivity
        std::map<std::string, double> sigma_params = {{"sigma_ref", 5.96e7}, {"alpha", 0.0039}, {"T_ref", 293.15}};
        copper.setTempDependentProperty("electrical_conductivity", sigma_params);

        copper.setProperty("thermal_conductivity", 400.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_2d_mesh(width, height, nx, ny)
            );
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Current2D>(copper));
        problem->addField(std::make_unique<Physics::Heat2D>(copper));
        problem->setup();
    }
};

// This test checks that internal Joule heating causes the temperature to rise over time.
TEST_F(TransientCoupled2DTest, InternalHeatingOverTime) {
    auto* heat_field = problem->getField("Temperature");
    auto* emag_field = problem->getField("Voltage");
    ASSERT_NE(heat_field, nullptr);
    ASSERT_NE(emag_field, nullptr);

    const double T_initial = 293.15;
    const double V_high = 0.01;

    // 1. Set Parameters
    problem->setTimeSteppingParameters(0.1, 1.0);

    // 2. Set Initial and Boundary Conditions
    heat_field->setInitialConditions(T_initial);

    auto& dof_manager = problem->getDofManager();
    const auto& mesh = problem->getMesh();
    for (const auto& node : mesh.getNodes()) {
        const auto& coords = node->getCoords();
        bool is_boundary = (std::abs(coords[0] - 0.0) < 1e-9 || std::abs(coords[0] - width) < 1e-9 ||
                            std::abs(coords[1] - 0.0) < 1e-9 || std::abs(coords[1] - height) < 1e-9);

        // Keep all boundaries at the initial temperature
        if (is_boundary) {
            heat_field->addBC(std::make_unique<Core::NeumannBC>
                (dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(0.1)));
        }

        // Apply voltage difference
        if (std::abs(coords[0] - 0.0) < 1e-9) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>
                (dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_high)));
        }
        else if (std::abs(coords[0] - width) < 1e-9) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>
                (dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(0.0)));
        }
    }

    // 3. Solve
    // problem->solveTransient();

    auto solver = Solver::SolverFactory::createSolver(*problem);
    solver->solveTransient(*problem);
    problem->exportResults("results_2d_transient_coupled.vtk");

    // 4. Validate
    double max_temp = heat_field->getSolution().maxCoeff();

    SimpleLogger::Logger::instance().info("Initial temperature: ", T_initial, " K");
    SimpleLogger::Logger::instance().info("Final max temperature: ", max_temp, " K");

    // The maximum temperature in the domain must be higher than the initial/boundary
    // temperature due to the internal heat generation.
    ASSERT_GT(max_temp, T_initial);
}
