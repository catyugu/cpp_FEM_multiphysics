#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>
#include <limits>
#include <map>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "physics/Current3D.hpp"
#include "physics/Heat3D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/coupling/ElectroThermalCoupling.hpp"

#undef max

class CoupledTransient3DTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    Core::Material copper{"Copper"};

    void SetUp() override {
        auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(0.02, 0.01, 0.01, 5, 5, 5));
        ASSERT_NE(mesh, nullptr);

        // Making conductivity dependent on temperature to properly test coupling
        copper.setProperty("electrical_conductivity", [](const std::map<std::string, double>& field_values) {
            double T = field_values.count("Temperature") ? field_values.at("Temperature") : 293.15;
            return 5.96e7 * (1.0 - 0.0039 * (T - 293.15)); // Simple linear model for copper resistivity
        });
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Current3D>(copper));
        problem->addField(std::make_unique<Physics::Heat3D>(copper));
        problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());

        problem->setTimeSteppingParameters(0.1, 1.0);
        problem->setIterativeSolverParameters(30, 1e-5);

        problem->setup();
    }
};

TEST_F(CoupledTransient3DTest, Simple3DTransientRun) {
    auto *emag_field = problem->getField("Voltage");
    auto *heat_field = problem->getField("Temperature");
    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    constexpr double V_in = 0.1;
    constexpr double T_initial = 293.15;
    constexpr double T_sink = 293.15;
    constexpr double bar_length_x = 0.02;
    constexpr double eps = 1e-8;

    heat_field->setInitialConditions(T_initial);

    // --- **FIX**: Define realistic boundary conditions ---
    auto bc_predicate = [&](const std::vector<double>& coords) {
        // Apply BCs only to the faces at x=0 and x=L
        return (std::abs(coords[0] - 0.0) < eps || std::abs(coords[0] - bar_length_x) < eps);
    };

    auto voltage_bc_value = [&](const std::vector<double>& coords) {
        return (std::abs(coords[0] - 0.0) < eps) ? V_in : 0.0;
    };

    auto heat_bc_value = [&](const std::vector<double>& /*coords*/) {
        return T_sink;
    };

    // Use the same predicate for both physics, applying BCs only to the ends
    auto voltage_bcs = Core::DirichletBC::create(problem->getDofManager(), problem->getMesh(), "Voltage", emag_field->getElementOrder(), bc_predicate, voltage_bc_value);
    auto heat_bcs = Core::DirichletBC::create(problem->getDofManager(), problem->getMesh(), "Temperature", heat_field->getElementOrder(), bc_predicate, heat_bc_value);

    emag_field->addBCs(std::move(voltage_bcs));
    heat_field->addBCs(std::move(heat_bcs));

    // Solve the transient problem
    ASSERT_NO_THROW(problem->solveTransient());
    problem->exportResults("results_coupled_transient_3D.vtk");

    // Validate
    const auto& final_temp_solution = heat_field->getSolution();
    double max_temp = final_temp_solution.maxCoeff();

    Utils::Logger::instance().info("Maximum final temperature: ", max_temp, " K");

    // The maximum temperature should now be noticeably higher than the sink temperature.
    ASSERT_GT(max_temp, T_sink + 1.0);
}