#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <fstream>
#include <cmath>
#include <map>
#include <limits>
#include <set> // Required for std::set

#include "core/Problem.hpp"
#include "core/Material.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include <solver/LinearSolver.hpp>
#include <utils/Exceptions.hpp>

#include "physics/Current3D.hpp"
#include "physics/Heat3D.hpp"
#include "io/Importer.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/coupling/ElectroThermalCoupling.hpp"
#include "core/mesh/TetElement.hpp" // For dynamic_cast to TetElement

#undef max

// Test fixture for validating 3D coupled simulation against a reference result file
class Coupled3DValidationTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    const std::string vtu_filename = "../data/electroThermalResults_3D.vtu";
    const std::string mesh_filename = "../data/electroThermalMesh_3D.mphtxt";
    Core::Material copper{"Copper"};

    void SetUp() override {
        // Setup problem
        std::unique_ptr<Core::Mesh> mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);
        ASSERT_NE(mesh, nullptr);

        copper.setProperty("electrical_conductivity", [](const std::map<std::string, double>& field_values) {
            double T = field_values.count("Temperature") ? field_values.at("Temperature") : 293.15;
            return 5.96e7 * (1.0 - 0.0039 * (T - 293.15)); // Simple linear model for copper resistivity
        });
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("thermal_capacity", 385.0);

        problem = std::make_unique<Core::Problem>(std::move(mesh));
        auto current_field = std::make_unique<Physics::Current3D>(copper);
        auto heat_field = std::make_unique<Physics::Heat3D>(copper);

        // Set element order to 2 for both fields
        // current_field->setElementOrder(2);
        // heat_field->setElementOrder(2);

        problem->addField(std::move(current_field));
        problem->addField(std::move(heat_field));
        problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());

        problem->setup();
        // problem->setLinearSolverType(Solver::SolverType::BiCGSTAB);
        problem->setIterativeSolverParameters(1000, 1e-3);
    }
};

TEST_F(Coupled3DValidationTest, CompareAgainstVtuResult) {
    constexpr double V_in = 1;
    constexpr double T_sink = 293.15;
    constexpr double bar_length = 1.0;
    constexpr double bar_width = 0.1;
    constexpr double bar_height = 0.1;
    constexpr double eps = 1e-5;

    auto *emag_field = problem->getField("Voltage");
    auto *heat_field = problem->getField("Temperature");

    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    // --- Define BCs using predicates and the new factory ---
    auto voltage_bc_predicate = [&](const std::vector<double>& coords) {
        return (std::abs(coords[0] - 0.0) < eps || std::abs(coords[0] - bar_length) < eps);
    };
    auto voltage_bc_value = [&](const std::vector<double>& coords) {
        return (std::abs(coords[0] - 0.0) < eps) ? V_in : 0.0;
    };
    auto heat_bc_predicate = [&](const std::vector<double>& coords) {
        return (std::abs(coords[0] - 0.0) < eps || std::abs(coords[0] - bar_length) < eps ||
                std::abs(coords[1] - 0.0) < eps || std::abs(coords[1] - bar_width) < eps ||
                std::abs(coords[2] - 0.0) < eps || std::abs(coords[2] - bar_height) < eps);
    };
    auto heat_bc_value = [&](const std::vector<double>& /*coords*/) {
        return T_sink;
    };

    auto voltage_bcs = Core::DirichletBC::create(problem->getDofManager(), problem->getMesh(), "Voltage", emag_field->getElementOrder(), voltage_bc_predicate, voltage_bc_value);
    auto heat_bcs = Core::DirichletBC::create(problem->getDofManager(), problem->getMesh(), "Temperature", heat_field->getElementOrder(), heat_bc_predicate, heat_bc_value);

    emag_field->addBCs(std::move(voltage_bcs));
    heat_field->addBCs(std::move(heat_bcs));


    // --- Solve ---
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results_3d.vtk");

    // --- Load Reference Data (Coordinates and Values) ---
    auto& logger = Utils::Logger::instance();
    // *** THIS IS THE FIX: Using the correct encoded names from the COMSOL VTU file ***
    std::vector<std::string> data_names = {"&#x6e29;&#x5ea6;", "&#x7535;&#x52bf;"}; // Names in COMSOL VTU
    IO::VtuData reference_data;
    try {
        reference_data = IO::Importer::read_vtu_points_and_data(vtu_filename, data_names);
    } catch (const Exception::FileIOException& e) {
        logger.error("Could not open reference VTU file: {}", vtu_filename);
        FAIL() << "Reference VTU file not found.";
    }

    // --- Robust Comparison ---
    const auto& fem_temp_solution = heat_field->getSolution();
    const auto& fem_volt_solution = emag_field->getSolution();

    ASSERT_TRUE(reference_data.point_data.count("&#x6e29;&#x5ea6;"));
    ASSERT_TRUE(reference_data.point_data.count("&#x7535;&#x52bf;"));

    const auto& ref_temp_values = reference_data.point_data.at("&#x6e29;&#x5ea6;");
    const auto& ref_volt_values = reference_data.point_data.at("&#x7535;&#x52bf;");

    double max_temp_diff = 0.0;
    double max_volt_diff = 0.0;
    int matched_nodes = 0;

    for (const auto& sim_node : problem->getMesh().getNodes()) {
        const auto& sim_coords = sim_node->getCoords();
        double min_dist_sq = std::numeric_limits<double>::max();
        size_t closest_ref_idx = -1;

        for (size_t i = 0; i < reference_data.points.size(); ++i) {
            double dist_sq = std::pow(sim_coords[0] - reference_data.points[i][0], 2) +
                             std::pow(sim_coords[1] - reference_data.points[i][1], 2) +
                             std::pow(sim_coords[2] - reference_data.points[i][2], 2);
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                closest_ref_idx = i;
            }
        }

        if (min_dist_sq < 1e-12) {
            matched_nodes++;
            int temp_dof = problem->getDofManager().getEquationIndex(sim_node->getId(), "Temperature");
            if (temp_dof != -1) {
                max_temp_diff = std::max(max_temp_diff, std::abs(fem_temp_solution(temp_dof) - ref_temp_values[closest_ref_idx]));
            }

            int volt_dof = problem->getDofManager().getEquationIndex(sim_node->getId(), "Voltage");
            if (volt_dof != -1) {
                max_volt_diff = std::max(max_volt_diff, std::abs(fem_volt_solution(volt_dof) - ref_volt_values[closest_ref_idx]));
            }
        }
    }

    logger.info("Successfully matched and compared ", matched_nodes," nodes based on coordinates.");
    logger.info("Maximum temperature difference: ", max_temp_diff, " K");
    logger.info("Maximum voltage difference: ",max_volt_diff," V");

    // Assert that the maximum differences are within an acceptable tolerance
    ASSERT_LT(max_temp_diff, 1);
    ASSERT_LT(max_volt_diff, 1e-5);
}