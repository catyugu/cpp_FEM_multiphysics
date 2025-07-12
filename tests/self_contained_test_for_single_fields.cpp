#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "physics/Magnetic1D.hpp"
#include "utils/SimpleLogger.hpp"

// Test fixture for 1D Magnetostatics
class Magnetic1DTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    Core::Material air{"Air"};

    void SetUp() override {
        // Define material properties for air
        air.setProperty("magnetic_permeability", 4.0 * EIGEN_PI * 1e-7); // mu_0

        // Create a simple 1D mesh
        auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_1d_mesh(1.0, 10));
        problem = std::make_unique<Core::Problem>(std::move(mesh));
        
        problem->addField(std::make_unique<Physics::Magnetic1D>(air));
        problem->setup();
    }
};

TEST_F(Magnetic1DTest, LinearPotential) {
    auto* magnetic_field = problem->getField("MagneticPotential");
    ASSERT_NE(magnetic_field, nullptr);

    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();
    
    // Apply Dirichlet boundary conditions to create a linear potential drop
    magnetic_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "MagneticPotential", Eigen::Vector<double, 1>(0.0)));
    magnetic_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 10, "MagneticPotential", Eigen::Vector<double, 1>(1.0)));

    // Solve the problem
    ASSERT_NO_THROW(problem->solveSteadyState());

    // Export results for visualization
    problem->exportResults("results_magnetic_1D.vtk");

    // Validate the solution
    const auto& solution = magnetic_field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        double x = node->getCoords()[0];
        double analytical_A = x; // Linear potential from 0 to 1
        double fem_A = solution(dof_manager.getEquationIndex(node->getId(), "MagneticPotential"));
        ASSERT_NEAR(fem_A, analytical_A, 1e-9);
    }
}