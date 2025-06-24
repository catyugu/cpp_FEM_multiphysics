#include "utils/SimpleLogger.hpp"
#include "core/Mesh.hpp"
#include "physics/EMag1D.hpp"
#include "physics/Heat1D.hpp"

#include <iostream>
#include <memory>

void solve_problem() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("--- Starting 1D Coupled FEM Solver ---");

    // 1. Define problem parameters
    const double length = 1.0; // m
    const int num_elements = 10;
    const double sigma = 5.96e7; // Electrical conductivity of copper
    const double k = 401.0;      // Thermal conductivity of copper

    // 2. Create the mesh
    std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));
    if(!mesh) {
        logger.error("Failed to create mesh.");
        return;
    }

    // 3. Setup the physics modules
    Physics::EMag1D emag_field(sigma);
    emag_field.setup(*mesh);

    Physics::Heat1D heat_field(k);
    heat_field.setup(*mesh);


    // 4. Assemble the systems (currently just logs the process)
    logger.info("\n--- Solving Electromagnetic Field ---");
    emag_field.assemble();
    // In the full solver, we would now:
    // a. Apply boundary conditions (e.g., V(0)=1V, V(L)=0V)
    // b. Solve the linear system for nodal voltages
    // c. Compute the electric field and Joule heat source in each element

    logger.info("\n--- Solving Heat Transfer Field ---");
    // d. Pass the computed heat source to the heat field
    //    std::vector<double> joule_heat_source = ...;
    //    heat_field.setHeatSource(joule_heat_source);
    heat_field.assemble();
    // e. Apply boundary conditions (e.g., T(0)=300K, T(L)=300K)
    // f. Solve the linear system for nodal temperatures

    logger.info("\n--- Simulation Finished ---");
}


int main() {
    // Configure logger
    SimpleLogger::Logger::instance().set_logfile("femsolver.log");
    SimpleLogger::Logger::instance().set_loglevel(SimpleLogger::LogLevel::info);

    try {
        solve_problem();
    } catch (const std::exception& e) {
        SimpleLogger::Logger::instance().error("An exception occurred: ", e.what());
        return 1;
    }

    return 0;
}
