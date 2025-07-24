#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Magnetic3D.hpp"
#include "core/sources/PrescribedCurrentDensity.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "utils/SimpleLogger.hpp"

// Test for 3D Magnetostatics
TEST(Magnetostatics3DTest, SolenoidTest) {
    auto& logger = Utils::Logger::instance();
    logger.info("--- Setting up 3D Magnetostatics Test: Solenoid ---");

    // 1. Create a simple mesh
    auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(1.0, 1.0, 5.0, 5, 5, 25));
    ASSERT_NE(mesh, nullptr);

    // 2. Define material
    Core::Material air("Air");
    air.setProperty("magnetic_permeability", 1.25663706e-6);

    // 3. Create problem and add the Magnetostatics3D field
    auto problem = std::make_unique<Core::Problem>(std::move(mesh));
    auto magnetic_field = std::make_unique<Physics::Magnetic3D>(air);
    // magnetic_field->setElementOrder(2);
    problem->addField(std::move(magnetic_field));

    problem->setup();

    auto* field = problem->getField("MagneticVectorPotential");
    ASSERT_NE(field, nullptr);

    // 4. Apply a circular current density to simulate a coil
    // This is a simplified approach using a uniform current in a region
    const double current_density_magnitude = 1.0e6; // 1 MA/m^2
    for (const auto& elem : problem->getMesh().getElements()) {
        // A simple way to create a "coil" region
        const auto& node_coords = elem->getNodes()[0]->getCoords();
        double r = std::sqrt(node_coords[0]*node_coords[0] + node_coords[1]*node_coords[1]);
        if (r > 0.4 && r < 0.5) {
            // Approximating circular current
            Eigen::Vector3d J(-node_coords[1]/r * current_density_magnitude, node_coords[0]/r * current_density_magnitude, 0);
            field->addSource(std::make_unique<Core::PrescribedCurrentDensity>(elem->getId(), J));
        }
    }

    // 5. Set A=0 on the outer boundary
    auto& dof_manager = problem->getDofManager();
    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        // Check if the node is on any of the 6 outer faces of the cube
        if (std::abs(coords[0] - 0.0) < 1e-9 || std::abs(coords[0] - 1.0) < 1e-9 ||
            std::abs(coords[1] - 0.0) < 1e-9 || std::abs(coords[1] - 1.0) < 1e-9 ||
            std::abs(coords[2] - 0.0) < 1e-9 || std::abs(coords[2] - 5.0) < 1e-9)
        {
            int base_dof_idx = dof_manager.getEquationIndex(node->getId(), "MagneticVectorPotential");
            if (base_dof_idx != -1) {
                // Create a separate BC for each component using the correct equation index
                field->addBC(std::make_unique<Core::DirichletBC>(base_dof_idx + 0, Eigen::Vector<double, 1>(0.0)));
                field->addBC(std::make_unique<Core::DirichletBC>(base_dof_idx + 1, Eigen::Vector<double, 1>(0.0)));
                field->addBC(std::make_unique<Core::DirichletBC>(base_dof_idx + 2, Eigen::Vector<double, 1>(0.0)));
            }
        }
    }
    // 6. Solve the problem
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("solenoid_results.vtk");

    // 7. Validate the results (Placeholder for B = curl(A) calculation)
    // Once a post-processor for calculating the curl of a vector field is available,
    // we will add assertions here to check for a uniform magnetic field inside the solenoid.
    SUCCEED() << "Solver ran without errors. Manual validation of solenoid_results.vtk is required.";
}