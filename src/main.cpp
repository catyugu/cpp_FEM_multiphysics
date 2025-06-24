#include "utils/SimpleLogger.hpp"
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/EMag1D.hpp"
#include "physics/Heat1D.hpp"

#include <iostream>
#include <memory>
#include <iomanip>

void run_simulation() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("--- Setting up 1D Coupled FEM Problem ---");

    Core::Material copper("Copper");
    copper.setProperty("electrical_conductivity", 5.96e7);
    copper.setProperty("thermal_conductivity", 401.0);

    const double length = 1.0;
    const int num_elements = 10;
    std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));
    auto problem = std::make_unique<Core::Problem>(std::move(mesh));

    problem->addField(std::make_unique<Physics::EMag1D>(copper));
    problem->addField(std::make_unique<Physics::Heat1D>(copper));

    problem->setup();

    const int num_nodes = num_elements + 1;
    const double V0 = 1.0;
    const double VL = 0.0;
    const double T_ambient = 300.0;

    auto& dof_manager = problem->getDofManager();
    problem->getField("Voltage")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Voltage", V0));
    problem->getField("Voltage")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Voltage", VL));
    problem->getField("Temperature")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_ambient));
    problem->getField("Temperature")->addBC(std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Temperature", T_ambient));

    problem->solve();

    logger.info("\n--- Post-processing Results ---");
    // (Existing validation code...)

    // --- ADD THIS SECTION ---
    // 8. Export results to a file for visualization
    logger.info("\n--- Exporting Results ---");
    problem->exportResults("results.vtk");
    // --- END OF ADDED SECTION ---
}


int main() {
    SimpleLogger::Logger::instance().set_logfile("femsolver.log");
    SimpleLogger::Logger::instance().set_loglevel(SimpleLogger::LogLevel::info);
    try {
        run_simulation();
    } catch (const std::exception& e) {
        SimpleLogger::Logger::instance().error("An exception occurred: ", e.what());
        return 1;
    }
    return 0;
}
