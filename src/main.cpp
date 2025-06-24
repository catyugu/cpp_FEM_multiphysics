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

    // 1. Define materials
    Core::Material copper("Copper");
    copper.setProperty("electrical_conductivity", 5.96e7);
    copper.setProperty("thermal_conductivity", 401.0);

    // 2. Define problem geometry and create problem
    const double length = 1.0;
    const int num_elements = 10;

    // FIX: Create the mesh first and wrap it in a unique_ptr
    std::unique_ptr<Core::Mesh> mesh(Core::Mesh::create_uniform_1d_mesh(length, num_elements));

    // Then, move ownership of the mesh to the Problem's constructor
    auto problem = std::make_unique<Core::Problem>(std::move(mesh));

    // 3. Add physics fields to the problem
    problem->addField(std::make_unique<Physics::EMag1D>(copper));
    problem->addField(std::make_unique<Physics::Heat1D>(copper));

    // 4. Setup the problem (builds DOF maps, sets up fields)
    problem->setup();

    // 5. Define and add boundary conditions
    const int num_nodes = num_elements + 1;
    const double V0 = 1.0;
    const double VL = 0.0;
    const double T_ambient = 300.0;

    auto& dof_manager = problem->getDofManager();
    problem->getField("Voltage")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, 0, "Voltage", V0)
    );
    problem->getField("Voltage")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Voltage", VL)
    );
    problem->getField("Temperature")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, 0, "Temperature", T_ambient)
    );
    problem->getField("Temperature")->addBC(
        std::make_unique<Core::DirichletBC>(dof_manager, num_nodes - 1, "Temperature", T_ambient)
    );

    // 6. Solve the problem
    problem->solve();

    // 7. Post-processing and validation
    logger.info("\n--- Post-processing Results ---");
    auto* heat_solution = &problem->getField("Temperature")->getSolution();

    logger.info("Heat Solution (Temperature):");
    std::cout << std::fixed << std::setprecision(4);
    for(int i = 0; i < num_nodes; ++i) {
        int eq_idx = dof_manager.getEquationIndex(i, "Temperature");
        logger.info("  Node ", i, ": ", (*heat_solution)(eq_idx), " K");
    }

    double max_temp_analytical = T_ambient + (copper.getProperty("electrical_conductivity") / (8.0 * copper.getProperty("thermal_conductivity"))) * std::pow(V0 - VL, 2);
    int center_node_idx = num_nodes / 2;
    int center_eq_idx = dof_manager.getEquationIndex(center_node_idx, "Temperature");
    double max_temp_fem = (*heat_solution)(center_eq_idx);

    logger.info("\n--- Validation ---");
    logger.info("Analytical max temperature: ", max_temp_analytical, " K");
    logger.info("FEM max temperature (center node): ", max_temp_fem, " K");
    logger.info("Difference: ", std::abs(max_temp_analytical - max_temp_fem), " K");
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
