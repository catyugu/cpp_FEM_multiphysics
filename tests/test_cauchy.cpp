#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"

class Cauchy1DTest : public ::testing::Test {
protected:
    const double length = 2.0;
    const int num_elements = 40;
    const int num_nodes = num_elements + 1;
    Core::Material aluminum{"Aluminum"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        aluminum.setProperty("thermal_conductivity", 237.0);
        aluminum.setProperty("density", 2700.0);
        aluminum.setProperty("specific_heat", 900.0);

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(
            length, num_elements));
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Heat1D>(aluminum));
        problem->setup();
    }
};

TEST_F(Cauchy1DTest, HeatConductionWithConvection) {
    const double T_fixed = 400.0;
    const double T_ambient = 293.15;
    const double h_conv = 15.0;

    auto& dof_manager = problem->getDofManager();
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", Eigen::Vector<double, 1>(T_fixed)));
    heat_field->addBC(std::make_unique<Core::CauchyBC>(dof_manager, num_nodes - 1, "Temperature",
        Eigen::Vector<double, 1>(h_conv), Eigen::Vector<double, 1>(T_ambient)));

    ASSERT_NO_THROW(problem->solveSteadyState());

    double k = aluminum.getProperty("thermal_conductivity");
    double C1 = (h_conv * (T_ambient - T_fixed)) / (h_conv * length + k);
    double C2 = T_fixed;

    const auto& solution_vector = heat_field->getSolution();

    for (int i = 0; i < num_nodes; ++i) {
        double x = problem->getMesh().getNode(i)->getCoords()[0];
        double analytical_T = C1 * x + C2;
        double fem_T = solution_vector(dof_manager.getEquationIndex(i, "Temperature"));
        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
}
