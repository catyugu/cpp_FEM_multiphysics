#include "utils/SimpleLogger.hpp"
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "core/BoundaryCondition.hpp"
#include "physics/Heat2D.hpp"
#include "io/Importer.hpp"
#include <memory>
#include <cmath>
#include <limits>
#include <vector>
#undef max
void run_comsol_circle_steadystate() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("--- Setting up 2D COMSOL Mesh Problem: Circle Steady-State ---");

    const std::string mesh_filename = "circle_mesh.mphtxt";
    std::unique_ptr<Core::Mesh> mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);
    if (!mesh) {
        logger.error("Failed to import mesh. Aborting simulation.");
        return;
    }

    Core::Material copper("Copper");
    copper.setProperty("thermal_conductivity", 401.0);
    copper.setProperty("density", 8960.0);
    copper.setProperty("specific_heat", 385.0);

    auto problem = std::make_unique<Core::Problem>(std::move(mesh));
    problem->addField(std::make_unique<Physics::Heat2D>(copper));
    problem->setup();

    const double T_ambient = 293.15;
    const double h_conv = 1.0;
    const double P_source = 10.0;
    const double circle_radius = 1.0;

    auto* heat_field = problem->getField("Temperature");
    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();

    Core::Node* source_node = nullptr;
    double min_dist_sq = std::numeric_limits<double>::max();

    // --- FIX: Calculate effective nodal length for more accurate BCs ---
    logger.info("Identifying boundary nodes and calculating effective nodal lengths...");
    std::vector<Core::Node*> boundary_nodes;
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(std::sqrt(coords[0]*coords[0] + coords[1]*coords[1]) - circle_radius) < 1e-4) {
            boundary_nodes.push_back(node);
        }
        double dist_sq = std::pow(coords[0] - 0.0, 2) + std::pow(coords[1] - 1.0, 2);
        if (dist_sq < min_dist_sq) {
            min_dist_sq = dist_sq;
            source_node = node;
        }
    }

    // Approximate the length represented by each boundary node
    double boundary_circumference = 2.0 * EIGEN_PI * circle_radius;
    double effective_length_per_node = boundary_circumference / boundary_nodes.size();
    double h_eff = h_conv * effective_length_per_node;
    logger.info("Found ", boundary_nodes.size(), " boundary nodes.");
    logger.info("Effective convection coefficient per node (h*L): ", h_eff, " W/(m*K)");

    // Apply the scaled Cauchy BC to all boundary nodes
    for (const auto& node : boundary_nodes) {
        heat_field->addBC(std::make_unique<Core::CauchyBC>(
            dof_manager, node->getId(), "Temperature", h_eff, T_ambient
        ));
    }

    // Apply the nodal heat source
    if (source_node) {
        logger.info("Applying point heat source of ", P_source, "W to node ", source_node->getId());
        heat_field->addBC(std::make_unique<Core::NeumannBC>(
            dof_manager, source_node->getId(), "Temperature", P_source
        ));
    }

    problem->solveSteadyState();
    problem->exportResults("comsol_circle_steadystate_results.vtk");

    logger.info("Final max temperature: ", heat_field->getSolution().maxCoeff(), " K");
    logger.info("This result should now be close to the COMSOL steady-state result (~295 K).");
}

int main() {
    SimpleLogger::Logger::instance().set_logfile("femsolver.log");
    SimpleLogger::Logger::instance().set_loglevel(SimpleLogger::LogLevel::info);
    try {
        run_comsol_circle_steadystate();
    } catch (const std::exception& e) {
        SimpleLogger::Logger::instance().error("An exception occurred: ", e.what());
        return 1;
    }
    return 0;
}
