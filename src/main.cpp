#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <core/bcs/BoundaryCondition.hpp>
#include "core/material/Material.hpp"
#include "core/Problem.hpp"
#include "io/Importer.hpp"
#include "physics/Heat2D.hpp"
#include "utils/Exceptions.hpp"
#include "utils/SimpleLogger.hpp"
#include "post/HeatFluxCalculator.hpp" // Include the new calculator

#undef max // Avoid conflict with std::numeric_limits<>::max() on Windows

void run_comsol_circle_steadystate() {
    auto& logger = Utils::Logger::instance();
    logger.info("--- Setting up 2D COMSOL Mesh Problem: Circle Steady-State ---");

    const std::string mesh_filename = "../data/circle_mesh.mphtxt";
    std::unique_ptr<Core::Mesh> mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);

    Core::Material copper("Copper");
    copper.setProperty("thermal_conductivity", 401.0);
    copper.setProperty("density", 8960.0);
    copper.setProperty("thermal_capacity", 385.0);

    auto problem = std::make_unique<Core::Problem>(std::move(mesh));
    problem->addField(std::make_unique<Physics::Heat2D>(copper));

    problem->getField("Temperature")->setElementOrder(2);

    // Add the new post-processor to the problem
    problem->addPostProcessor(std::make_unique<Post::HeatFluxCalculator>());

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

    double boundary_circumference = 2.0 * 3.14159265358979323846 * circle_radius;
    double effective_length_per_node = boundary_circumference / boundary_nodes.size();
    double h_eff = h_conv * effective_length_per_node;
    logger.info("Found ", boundary_nodes.size(), " boundary nodes.");
    logger.info("Effective convection coefficient per node (h*L): ", h_eff, " W/(m*K)");

    for (const auto& node : boundary_nodes) {
        heat_field->addBC(std::make_unique<Core::CauchyBC>(
            dof_manager, node->getId(), "Temperature", Eigen::Vector<double, 1>(h_eff), Eigen::Vector<double, 1>(T_ambient)
        ));
    }

    if (source_node) {
        logger.info("Applying point heat source of ", P_source, "W to node ", source_node->getId());
        heat_field->addBC(std::make_unique<Core::NeumannBC>(
            dof_manager, source_node->getId(), "Temperature", Eigen::Vector<double, 1>(P_source)
        ));
    }

    problem->solveSteadyState();
    problem->exportResults("comsol_circle_steadystate_results.vtk");

    logger.info("Final max temperature: ", heat_field->getSolution().maxCoeff(), " K");
}

int main() {
    Utils::Logger::instance().set_logfile("femsolver.log");
    Utils::Logger::instance().set_loglevel(Utils::LogLevel::info);
    try {
        run_comsol_circle_steadystate();
    } catch (const Exception::FileIOException& e) {
        Utils::Logger::instance().error("A file I/O error occurred: ", e.what());
        return 1;
    } catch (const Exception::SolverException& e) {
        Utils::Logger::instance().error("A solver error occurred: ", e.what());
        return 1;
    } catch (const Exception::ConfigurationException& e) {
        Utils::Logger::instance().error("A configuration error occurred: ", e.what());
        return 1;
    }
    catch (const std::exception& e) {
        Utils::Logger::instance().error("An unexpected error occurred: ", e.what());
        return 1;
    }
    return 0;
}
