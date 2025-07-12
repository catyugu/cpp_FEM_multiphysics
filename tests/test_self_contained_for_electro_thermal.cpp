#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <cmath>
#include <limits>

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "core/ElectroThermalCoupling.hpp"
#include "physics/Heat1D.hpp"
#include "physics/Heat2D.hpp"
#include "physics/Heat3D.hpp"
#include "physics/Current1D.hpp"
#include "physics/Current2D.hpp"
#include "io/Importer.hpp"
#include "utils/SimpleLogger.hpp"

// =================================================================================================
// Test Fixture for 1D Single Field Heat Transfer
// =================================================================================================
class Heat1DTest : public ::testing::Test {
protected:
    const double length = 2.0;
    const int num_elements = 40;
    const int num_nodes = num_elements + 1;
    Core::Material material{"TestMaterial"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        material.setProperty("thermal_conductivity", 50.0);
        material.setProperty("density", 7850.0);
        material.setProperty("specific_heat", 462.0);

        std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Heat1D>(material));
        problem->setup();
    }
};

TEST_F(Heat1DTest, NeumannBoundaryCondition) {
    const double T_fixed = 373.15;
    const double flux_out = 1000.0;

    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);
    auto& dof_manager = problem->getDofManager();

    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", Eigen::Vector<double, 1>(T_fixed)));
    heat_field->addBC(std::make_unique<Core::NeumannBC>(dof_manager, num_nodes - 1, "Temperature", -Eigen::Vector<double, 1>(flux_out)));

    ASSERT_NO_THROW(problem->solveSteadyState());

    double k = material.getProperty("thermal_conductivity");
    const auto& solution = heat_field->getSolution();

    for (int i = 0; i < num_nodes; ++i) {
        double x = problem->getMesh().getNode(i)->getCoords()[0];
        double analytical_T = T_fixed - (flux_out / k) * x;
        double fem_T = solution(dof_manager.getEquationIndex(i, "Temperature"));
        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
}

TEST_F(Heat1DTest, CauchyBoundaryCondition) {
    const double T_fixed = 400.0;
    const double T_ambient = 293.15;
    const double h_conv = 15.0;

    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);
    auto& dof_manager = problem->getDofManager();

    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", Eigen::Vector<double, 1>(T_fixed)));
    heat_field->addBC(std::make_unique<Core::CauchyBC>(dof_manager, num_nodes - 1, "Temperature", Eigen::Vector<double, 1>(h_conv), Eigen::Vector<double, 1>(T_ambient)));

    ASSERT_NO_THROW(problem->solveSteadyState());

    double k = material.getProperty("thermal_conductivity");
    double C1 = (h_conv * (T_ambient - T_fixed)) / (h_conv * length + k);
    double C2 = T_fixed;
    const auto& solution = heat_field->getSolution();

    for (int i = 0; i < num_nodes; ++i) {
        double x = problem->getMesh().getNode(i)->getCoords()[0];
        double analytical_T = C1 * x + C2;
        double fem_T = solution(dof_manager.getEquationIndex(i, "Temperature"));
        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
}

// =================================================================================================
// Test Fixture for 2D/3D Single Field and Matrix Validation
// =================================================================================================
class HeatMultiDTest : public ::testing::Test {
protected:
    Core::Material material{"TestMaterial"};
};

TEST_F(HeatMultiDTest, Heat2DPlateFixedBoundaries) {
    const double size = 1.0;
    const int num_div = 10;
    material.setProperty("thermal_conductivity", 401.0);
    material.setProperty("density", 8960.0);
    material.setProperty("specific_heat", 385.0);

    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(
        Core::Mesh::create_uniform_2d_mesh(size, size, num_div, num_div)
    ));
    problem->addField(std::make_unique<Physics::Heat2D>(material));
    problem->setup();

    const double T_hot = 500.0, T_cold = 300.0;
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[1] - size) < 1e-9) { // Top edge
            heat_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Temperature", Eigen::Vector<double, 1>(T_hot)));
        } else if (std::abs(coords[0] - 0.0) < 1e-9 || std::abs(coords[0] - size) < 1e-9 || std::abs(coords[1] - 0.0) < 1e-9) { // Other edges
            heat_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Temperature", Eigen::Vector<double, 1>(T_cold)));
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_2d_heat.vtk");

    const auto& solution = heat_field->getSolution();
    ASSERT_LE(solution.maxCoeff(), T_hot + 1e-9);
    ASSERT_GE(solution.minCoeff(), T_cold - 1e-9);
}

TEST_F(HeatMultiDTest, Heat3DCubeFixedBoundaries) {
    const double size = 1.0;
    const int num_div = 5;
    material.setProperty("thermal_conductivity", 50.2);
    material.setProperty("density", 7850.0);
    material.setProperty("specific_heat", 462.0);

    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(
        Core::Mesh::create_uniform_3d_mesh(size, size, size, num_div, num_div, num_div)
    ));
    problem->addField(std::make_unique<Physics::Heat3D>(material));
    problem->setup();

    const double T_hot = 400.0, T_cold = 300.0;
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) { // Left face
            heat_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Temperature", Eigen::Vector<double, 1>(T_hot)));
        } else if (std::abs(coords[0] - size) < 1e-9) { // Right face
            heat_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Temperature", Eigen::Vector<double, 1>(T_cold)));
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_3d_heat.vtk");

    const auto& solution = heat_field->getSolution();
    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        double analytical_T = T_hot - (T_hot - T_cold) * (coords[0] / size);
        double fem_T = solution(problem->getDofManager().getEquationIndex(node->getId(), "Temperature"));
        ASSERT_NEAR(fem_T, analytical_T, 1e-9);
    }
}

TEST_F(HeatMultiDTest, MassMatrixSummation) {
    const double width = 0.5, height = 0.2;
    material.setProperty("thermal_conductivity", 1.0);
    material.setProperty("density", 100.0);
    material.setProperty("specific_heat", 2.0);

    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(
        Core::Mesh::create_uniform_2d_mesh(width, height, 2, 1)
    ));
    problem->addField(std::make_unique<Physics::Heat2D>(material));
    problem->setup();

    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);
    heat_field->assemble();

    const auto& M = heat_field->getMassMatrix();
    double total_mass = M.sum();
    double expected_total_mass = material.getProperty("density") * material.getProperty("specific_heat") * (width * height);

    ASSERT_NEAR(total_mass, expected_total_mass, 1e-9);
}


// =================================================================================================
// Test Fixture for 2D Single Field Current
// =================================================================================================
class Current2DTest : public ::testing::Test {
protected:
    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;

    void SetUp() override {
        copper.setProperty("electrical_conductivity", 5.96e7);
        auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_2d_mesh(2.0, 1.0, 20, 10));
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        problem->addField(std::make_unique<Physics::Current2D>(copper));
        problem->setup();
    }
};

TEST_F(Current2DTest, RectangularPlateVoltageDrop) {
    const double V_high = 5.0, V_low = 1.0, width = 2.0;
    auto* emag_field = problem->getField("Voltage");
    ASSERT_NE(emag_field, nullptr);

    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Voltage", Eigen::Vector<double, 1>(V_high)));
        } else if (std::abs(coords[0] - width) < 1e-9) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Voltage", Eigen::Vector<double, 1>(V_low)));
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_2d_current.vtk");

    const auto& solution = emag_field->getSolution();
    for (const auto& node : problem->getMesh().getNodes()) {
        double x = node->getCoords()[0];
        double analytical_V = V_high - (V_high - V_low) * (x / width);
        double fem_V = solution(problem->getDofManager().getEquationIndex(node->getId(), "Voltage"));
        ASSERT_NEAR(fem_V, analytical_V, 1e-9);
    }
}

// =================================================================================================
// Test Fixture for Coupled Electro-Thermal Problems
// =================================================================================================
class CoupledTest : public ::testing::Test {
protected:
    Core::Material copper{"Copper"};
};

TEST_F(CoupledTest, Coupled1DSteadyState) {
    const double length = 1.0, V0 = 1.0, T_ambient = 300.0;
    copper.setProperty("electrical_conductivity", 5.96e7);
    copper.setProperty("thermal_conductivity", 401.0);
    copper.setProperty("density", 8960.0);
    copper.setProperty("specific_heat", 385.0);

    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(
        Core::Mesh::create_uniform_1d_mesh(length, 20)
    ));
    problem->addField(std::make_unique<Physics::Current1D>(copper));
    problem->addField(std::make_unique<Physics::Heat1D>(copper));
    problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());
    problem->setup();

    auto* emag_field = problem->getField("Voltage");
    auto* heat_field = problem->getField("Temperature");
    auto& dof_manager = problem->getDofManager();

    emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Voltage", Eigen::Vector<double, 1>(V0)));
    emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 20, "Voltage", Eigen::Vector<double, 1>(0.0)));
    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", Eigen::Vector<double, 1>(T_ambient)));
    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 20, "Temperature", Eigen::Vector<double, 1>(T_ambient)));

    ASSERT_NO_THROW(problem->solveSteadyState());

    double k = copper.getProperty("thermal_conductivity");
    double analytical_max_temp = T_ambient + (copper.getProperty("electrical_conductivity") / (8.0 * k)) * std::pow(V0, 2);
    double fem_max_temp = heat_field->getSolution().maxCoeff();

    ASSERT_NEAR(fem_max_temp, analytical_max_temp, analytical_max_temp * 0.01);
}

TEST_F(CoupledTest, Coupled2DTransient) {
    const double width = 0.02, height = 0.01, T_initial = 293.15, V_high = 0.01;
    std::map<std::string, double> sigma_params = {{"sigma_ref", 5.96e7}, {"alpha", 0.0039}, {"T_ref", 293.15}};
    copper.setTempDependentProperty("electrical_conductivity", sigma_params);
    copper.setProperty("thermal_conductivity", 400.0);
    copper.setProperty("density", 8960.0);
    copper.setProperty("specific_heat", 385.0);

    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(
        Core::Mesh::create_uniform_2d_mesh(width, height, 20, 10)
    ));
    problem->addField(std::make_unique<Physics::Current2D>(copper));
    problem->addField(std::make_unique<Physics::Heat2D>(copper));
    problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());
    problem->setup();

    problem->setTimeSteppingParameters(0.1, 1.0);
    auto* heat_field = problem->getField("Temperature");
    auto* emag_field = problem->getField("Voltage");
    heat_field->setInitialConditions(T_initial);

    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) { // Left edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Voltage", Eigen::Vector<double, 1>(V_high)));
            heat_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Temperature", Eigen::Vector<double, 1>(T_initial)));
        } else if (std::abs(coords[0] - width) < 1e-9) { // Right edge
            emag_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Voltage", Eigen::Vector<double, 1>(0.0)));
            heat_field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), node->getId(), "Temperature", Eigen::Vector<double, 1>(T_initial)));
        }
    }

    ASSERT_NO_THROW(problem->solveTransient());
    problem->exportResults("results_2d_transient_coupled.vtk");

    ASSERT_GT(heat_field->getSolution().maxCoeff(), T_initial);
}


// =================================================================================================
// Test Fixture for COMSOL Mesh Benchmarks
// =================================================================================================
class ComsolBenchmarkTest : public ::testing::Test {
protected:
    Core::Material copper{"Copper"};
    std::unique_ptr<Core::Problem> problem;
    const std::string mesh_filename = "busbar_mesh.mphtxt";

    void SetUp() override {
        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        auto mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);
        ASSERT_NE(mesh, nullptr);
        problem = std::make_unique<Core::Problem>(std::move(mesh));
    }
};

TEST_F(ComsolBenchmarkTest, CoupledBusbarVsComsol) {
    const double COMSOL_MAX_TEMP = 480.0; // K (Reference value)
    const double V_in = 0.1, T_sink = 293.15, bar_width = 0.1;

    problem->addField(std::make_unique<Physics::Current2D>(copper));
    problem->addField(std::make_unique<Physics::Heat2D>(copper));
    problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());
    problem->setup();

    auto* emag_field = problem->getField("Voltage");
    auto* heat_field = problem->getField("Temperature");
    auto& dof_manager = problem->getDofManager();

    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-4) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(V_in)));
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_sink)));
        } else if (std::abs(coords[0] - bar_width) < 1e-4) {
            emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(0.0)));
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_sink)));
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("busbar_results.vtk");

    double fem_max_temp = heat_field->getSolution().maxCoeff();
    SimpleLogger::Logger::instance().info("COMSOL Max Temperature: ", COMSOL_MAX_TEMP, " K");
    SimpleLogger::Logger::instance().info("Our FEM Max Temperature: ", fem_max_temp, " K");

    ASSERT_NEAR(fem_max_temp, COMSOL_MAX_TEMP, COMSOL_MAX_TEMP * 0.01);
}
