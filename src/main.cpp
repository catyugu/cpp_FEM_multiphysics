#include <iostream>
#include <memory>
#include <cmath>
#include "SimpleLogger.hpp"
#include "mesh.h"
#include "HeatField.h"
SimpleLogger::Logger& logger = SimpleLogger::Logger::instance();
// The main function to run the FEM heat simulation
int main() {
    logger.info("Starting FEM Heat Solver...");

    // 1. Load the mesh
    auto mesh_ptr = std::make_unique<Mesh>();
    if (!mesh_ptr->load_mesh("../my_mesh.txt")) {
        logger.error("Mesh loading failed. Exiting.");
        return -1;
    }
    // 2. Initialize the Heat Field with the loaded mesh
    HeatField heat_field(std::move(mesh_ptr));

    // 3. Define and Apply Simulation Conditions

    // Initial condition: entire body at 20 degrees C
    heat_field.apply_init_condition([](const Eigen::Vector3d& coord) {
        if (std::abs(coord.x()) < 0.001) {
            return 80.0; // Hot boundary
        }
        return 20.0; // Cold boundary
    });

    // Boundary condition: 100 degrees on the left edge (x=0), 20 degrees elsewhere
    // heat_field.apply_border_condition([](const Eigen::Vector3d& coord) {
    //     if (std::abs(coord.x()) < 1e-6) {
    //         return 100.0; // Hot boundary
    //     }
    //     return 20.0; // Cold boundary
    // });

    // No external heat source in this example
    // heat_field.apply_extern_field(...)

    // 4. Run the simulation
    const double time_step = 1;
    const int num_steps = 30;
    logger.info("Running simulation for " + std::to_string(num_steps) + " steps with dt=" + std::to_string(time_step));

    for (int i = 0; i < num_steps; ++i) {
        heat_field.step_forward(time_step);
        logger.info("Completed step " + std::to_string(i+1) + "/" + std::to_string(num_steps));
    }

    // 5. Write the final result to a VTK file
    heat_field.write_vtk("result.vtk");

    logger.info("Simulation finished. Result saved to result.vtk");

    return 0;
}