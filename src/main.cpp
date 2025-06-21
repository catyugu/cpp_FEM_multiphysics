#include "Analysis.h"
#include "mesh.h"
#include <memory>
#include <cmath>

// Define the global logger instance
SimpleLogger::Logger& logger = SimpleLogger::Logger::instance();

int main() {
    logger.info("Starting Refactored FEM Solver...");

    // 1. Load the mesh
    auto mesh_ptr = std::make_unique<Mesh>(4);
    if (!mesh_ptr->load_mesh("../my_mesh.txt")) {
        logger.error("Mesh loading failed. Exiting.");
        return -1;
    }

    // 2. Create the Analysis object
    Analysis simulation(std::move(mesh_ptr));

    // 3. Setup the simulation
    simulation.initialize_physics(); // Important: Creates the HeatTetrahedron objects

    simulation.apply_initial_condition([](const Eigen::Vector3d& coord) {
        if (std::abs(coord.x()) < 0.01) return 100.0;
        return 20.0; // Initial temperature of 20 C
    });

    // simulation.apply_boundary_condition([](const Eigen::Vector3d& coord) {
    //     if (std::abs(coord.x()) < 0.01) return 100.0; // Hot boundary at x=0
    //     return 20.0; // Cold boundary elsewhere
    // });

    // 4. Run the simulation and save results
    simulation.run_simulation(10, 0.1);
    simulation.write_vtk("result.vtk");

    logger.info("Simulation finished.");
    return 0;
}