#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "utils/SimpleLogger.hpp"
#include "physics/Current3D.hpp"
#include "physics/Magnetic3D.hpp"
#include "physics/Heat2D.hpp"
#include <cmath> // For std::sin, std::sinh
#include <limits> // For std::numeric_limits
#undef max

TEST(Current3DTest, UniformCube) {
    Core::Material material("TestMaterial");
    material.setProperty("electrical_conductivity", 1.0);

    auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(1.0, 1.0, 1.0, 2, 2, 2));
    auto problem = std::make_unique<Core::Problem>(std::move(mesh));
    problem->addField(std::make_unique<Physics::Current3D>(material));
    problem->setup();

    auto* current_field = problem->getField("Voltage");
    ASSERT_NE(current_field, nullptr);

    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();

    for(const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) {
            current_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(0.0)));
        }
        if (std::abs(coords[0] - 1.0) < 1e-9) {
            current_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Voltage", Eigen::Vector<double, 1>(1.0)));
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());

    const auto& solution = current_field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        double x = node->getCoords()[0];
        double analytical_V = x;
        double fem_V = solution(dof_manager.getEquationIndex(node->getId(), "Voltage"));
        ASSERT_NEAR(fem_V, analytical_V, 1e-9);
    }
}

TEST(Magnetic3DTest, UniformCube) {
    Core::Material material("TestMaterial");
    material.setProperty("magnetic_permeability", 1.0);

    auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(1.0, 1.0, 1.0, 2, 2, 2));
    auto problem = std::make_unique<Core::Problem>(std::move(mesh));
    problem->addField(std::make_unique<Physics::Magnetic3D>(material));
    problem->setup();

    auto* magnetic_field = problem->getField("MagneticPotential");
    ASSERT_NE(magnetic_field, nullptr);

    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();

    for(const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9) {
            magnetic_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "MagneticPotential", Eigen::Vector<double, 1>(0.0)));
        }
        if (std::abs(coords[0] - 1.0) < 1e-9) {
            magnetic_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "MagneticPotential", Eigen::Vector<double, 1>(1.0)));
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());

    const auto& solution = magnetic_field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        double x = node->getCoords()[0];
        double analytical_A = x;
        double fem_A = solution(dof_manager.getEquationIndex(node->getId(), "MagneticPotential"));
        ASSERT_NEAR(fem_A, analytical_A, 1e-9);
    }
}
class Heat2DValidationTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    Core::Material testMaterial{"TestMaterial"};
    const double tolerance = 1e-3;

    void SetUp() override {
        // Define material properties
        testMaterial.setProperty("thermal_conductivity", 1.0);
        testMaterial.setProperty("density", 1.0);
        testMaterial.setProperty("specific_heat", 1.0);

        // Create a simple 2D mesh
        auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_2d_mesh(1.0, 1.0, 20, 20));
        problem = std::make_unique<Core::Problem>(std::move(mesh));

        problem->addField(std::make_unique<Physics::Heat2D>(testMaterial));
        problem->setup();
    }
};

// Test against a known analytical solution for a square plate
TEST_F(Heat2DValidationTest, SquarePlateAnalyticalSolution) {
    auto* heat_field = problem->getField("Temperature");
    ASSERT_NE(heat_field, nullptr);

    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();
    const double pi = EIGEN_PI;

    // Apply Dirichlet boundary conditions
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        double x = coords[0];
        double y = coords[1];

        // T=0 on left, right, and bottom edges
        if (std::abs(x - 0.0) < 1e-9 || std::abs(x - 1.0) < 1e-9 || std::abs(y - 0.0) < 1e-9) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(0.0)));
        }
        // T=sin(pi*x) on top edge
        if (std::abs(y - 1.0) < 1e-9) {
            double T_val = std::sin(pi * x);
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_val)));
        }
    }

    // Solve the steady-state problem
    ASSERT_NO_THROW(problem->solveSteadyState());

    // Validate the solution against the analytical result
    const auto& solution = heat_field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        double x = coords[0];
        double y = coords[1];

        // Analytical solution: T(x,y) = sin(pi*x) * sinh(pi*y) / sinh(pi)
        double analytical_T = std::sin(pi * x) * std::sinh(pi * y) / std::sinh(pi);

        int dof_idx = dof_manager.getEquationIndex(node->getId(), "Temperature");
        double fem_T = solution(dof_idx);

        ASSERT_NEAR(fem_T, analytical_T, tolerance);
    }

    SimpleLogger::Logger::instance().info("Heat2D standalone test passed against analytical solution.");
}




// Test fixture for validating 2D Heat Transfer with higher-order quadrature
class Heat2DHigherOrderQuadratureTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    Core::Material testMaterial{"TestMaterial"};
    const double tolerance = 1e-3; // Tolerance for comparison against analytical solution

    void SetUp() override {
        // Define material properties
        testMaterial.setProperty("thermal_conductivity", 1.0); //
        testMaterial.setProperty("density", 1.0);
        testMaterial.setProperty("specific_heat", 1.0);

        // Create a simple 2D mesh (linear elements)
        auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_2d_mesh(1.0, 1.0, 20, 20)); //
        problem = std::make_unique<Core::Problem>(std::move(mesh)); //

        // Add the Heat2D physics field
        problem->addField(std::make_unique<Physics::Heat2D>(testMaterial)); //

        // Get the heat field and set its element order for quadrature
        auto* heat_field = problem->getField("Temperature"); //
        ASSERT_NE(heat_field, nullptr);
        heat_field->setElementOrder(5); // Set to a higher order for quadrature (e.g., 3)

        problem->setup(); //
    }
};

// Test against a known analytical solution for a square plate, using higher-order quadrature
TEST_F(Heat2DHigherOrderQuadratureTest, SquarePlateAnalyticalSolutionWithHighOrderQuadrature) {
    auto* heat_field = problem->getField("Temperature"); //
    ASSERT_NE(heat_field, nullptr);

    auto& dof_manager = problem->getDofManager(); //
    const auto& mesh_ref = problem->getMesh(); //

    // Apply Dirichlet boundary conditions
    for (const auto& node : mesh_ref.getNodes()) { //
        const auto& coords = node->getCoords(); //
        double x = coords[0];
        double y = coords[1];

        // T=0 on left (x=0), right (x=1), and bottom (y=0) edges
        if (std::abs(x - 0.0) < 1e-9 || std::abs(x - 1.0) < 1e-9 || std::abs(y - 0.0) < 1e-9) {
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(0.0))); //
        }
        // T=sin(pi*x) on top edge (y=1)
        if (std::abs(y - 1.0) < 1e-9) {
            double T_val = std::sin(EIGEN_PI * x);
            heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(T_val))); //
        }
    }

    // Solve the steady-state problem
    ASSERT_NO_THROW(problem->solveSteadyState()); //

    // Validate the solution against the analytical result
    const auto& solution = heat_field->getSolution(); //
    for (const auto& node : mesh_ref.getNodes()) { //
        const auto& coords = node->getCoords(); //
        double x = coords[0];
        double y = coords[1];

        // Analytical solution: T(x,y) = sin(pi*x) * sinh(pi*y) / sinh(pi)
        double analytical_T = std::sin(EIGEN_PI * x) * std::sinh(EIGEN_PI * y) / std::sinh(EIGEN_PI);

        int dof_idx = dof_manager.getEquationIndex(node->getId(), "Temperature"); //
        ASSERT_NE(dof_idx, -1); // Ensure DOF exists for this node
        double fem_T = solution(dof_idx);

        ASSERT_NEAR(fem_T, analytical_T, tolerance);
    }

    SimpleLogger::Logger::instance().info("Heat2D Higher-Order Quadrature test passed against analytical solution."); //
}