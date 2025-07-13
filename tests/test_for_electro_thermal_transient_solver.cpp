#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>
#include <limits>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "physics/Current3D.hpp" // Use Current3D
#include "physics/Heat3D.hpp"    // Use Heat3D
#include "utils/SimpleLogger.hpp"
#include "core/coupling/ElectroThermalCoupling.hpp"

#undef max // Undefine max macro potentially defined by Windows headers

class CoupledTransient3DTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    Core::Material copper{"Copper"};

    void SetUp() override {
        // Setup problem with a uniform 3D mesh
        // Using a smaller mesh for faster test execution, e.g., 2x2x2 hexes split into tets
        auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(0.02, 0.01, 0.01, 10, 10, 10)); // Smaller cube
        ASSERT_NE(mesh, nullptr);

        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Current3D>(copper)); // Use 3D physics
        problem->addField(std::make_unique<Physics::Heat3D>(copper));    // Use 3D physics
        problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());

        // Set time stepping parameters
        problem->setTimeSteppingParameters(0.1, 1.0); // 0.1s time step, 1.0s total time
        problem->setIterativeSolverParameters(20, 1e-4); // Max 20 inner iterations, 1e-4 tolerance

        problem->setup();
    }
};

TEST_F(CoupledTransient3DTest, Simple3DTransientRun) {
    auto *emag_field = problem->getField("Voltage");
    auto *heat_field = problem->getField("Temperature");
    auto &dof_manager = problem->getDofManager();
    const auto &mesh_ref = problem->getMesh();
    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    constexpr double V_in = 0.1;
    constexpr double T_initial = 293.15; // Initial temperature in K
    constexpr double T_sink = 293.15;
    constexpr double bar_length_x = 0.02;
    constexpr double bar_width_y = 0.01;
    constexpr double bar_height_z = 0.01;
    constexpr double eps = 1e-8;

    // Set initial conditions for temperature
    heat_field->setInitialConditions(T_initial);

    // Apply boundary conditions
    for (const auto &node: mesh_ref.getNodes()) {
        const auto &coords = node->getCoords();
        // Voltage BCs: 0.1V at x=0, 0V at x=bar_length_x
        if (std::abs(coords[0] - 0.0) < eps) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_in)));
        }
        if (std::abs(coords[0] - bar_length_x) < eps) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(0.0)));
        }
        // Heat BCs: Fixed temperature on all outer surfaces (x, y, z min/max)
        if (std::abs(coords[0] - 0.0) < eps || std::abs(coords[0] - bar_length_x) < eps ||
            std::abs(coords[1] - 0.0) < eps || std::abs(coords[1] - bar_width_y) < eps ||
            std::abs(coords[2] - 0.0) < eps || std::abs(coords[2] - bar_height_z) < eps) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_sink)));
        }
    }

    // Store initial solution for comparison
    Eigen::MatrixXd initial_temp_solution = heat_field->getSolution();

    // Solve the transient problem
    ASSERT_NO_THROW(problem->solveTransient());

    problem->exportResults("results_coupled_transient_3D.vtk");

    // Validate: Check if temperature has changed from initial state
    const auto& final_temp_solution = heat_field->getSolution();

    double max_temp_diff = 0.0;
    for (int i = 0; i < final_temp_solution.size(); ++i) {
        max_temp_diff = std::max(max_temp_diff, std::abs(final_temp_solution(i) - initial_temp_solution(i)));
    }

    SimpleLogger::Logger::instance().info("Maximum temperature change from initial: ", max_temp_diff, " K");

    // Expect some temperature change if Joule heating occurs and boundary conditions allow
    ASSERT_GT(max_temp_diff, 1e-6); // Expect a noticeable change
}