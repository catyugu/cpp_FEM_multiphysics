#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Current1D.hpp"
#include "physics/Heat1D.hpp"
#include "core/BoundaryCondition.hpp"
#include "utils/SimpleLogger.hpp"

class TransientCoupled1DTest : public ::testing::Test {
protected:
    const double length = 0.1;
    const int num_elements = 10;
    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        std::map<std::string, double> sigma_params = {{"sigma_ref", 5.96e7}, {"alpha", 0.0039}, {"T_ref", 293.15}};
        copper.setTempDependentProperty("electrical_conductivity", sigma_params);
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        std::unique_ptr<Core::Mesh> mesh(
            Core::Mesh::create_uniform_1d_mesh(length, num_elements));
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Current1D>(copper));
        problem->addField(std::make_unique<Physics::Heat1D>(copper));
        problem->setup();
    }
};

TEST_F(TransientCoupled1DTest, RodHeatingOverTime) {
    auto* heat_field = problem->getField("Temperature");
    auto* emag_field = problem->getField("Voltage");
    const double T_initial = 300.0, V_high = 0.5;

    problem->setTimeSteppingParameters(0.1, 1.0);
    heat_field->setInitialConditions(T_initial);

    auto& dof_manager = problem->getDofManager();
    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_initial));
    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_elements, "Temperature", T_initial));
    emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Voltage", V_high));
    emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_elements, "Voltage", 0.0));

    problem->solveTransient();
    problem->exportResults("results_1d_transient_coupled.vtk");

    double max_temp = heat_field->getSolution().maxCoeff();
    ASSERT_GT(max_temp, T_initial);
}
