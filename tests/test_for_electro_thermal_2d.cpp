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
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "io/Importer.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/coupling/ElectroThermalCoupling.hpp"
#include "core/mesh/TriElement.hpp" // For dynamic_cast to TriElement

#undef max

// Test fixture for validating coupled simulation against a reference result file
class Coupled2DValidationTest : public ::testing::Test {
protected:
    std::unique_ptr<Core::Problem> problem;
    const std::string mesh_filename = "../data/electroThermalMesh_2D.mphtxt";
    const std::string vtu_filename = "../data/electroThermalResults_2D.vtu";
    Core::Material copper{"Copper"};

    void SetUp() override {
        // Setup problem
        std::unique_ptr<Core::Mesh> mesh = IO::Importer::read_comsol_mphtxt(mesh_filename);
        ASSERT_NE(mesh, nullptr);

        copper.setProperty("electrical_conductivity", 5.96e7);
        copper.setProperty("thermal_conductivity", 401.0);
        copper.setProperty("density", 8960.0);
        copper.setProperty("specific_heat", 385.0);

        problem = std::make_unique<Core::Problem>(std::move(mesh));
        auto current_field = std::make_unique<Physics::Current2D>(copper);
        auto heat_field = std::make_unique<Physics::Heat2D>(copper);

        // Set element order to 2 for both fields
        current_field->setElementOrder(2);
        heat_field->setElementOrder(2);

        problem->addField(std::move(current_field));
        problem->addField(std::move(heat_field));
        problem->getCouplingManager().addCoupling(std::make_unique<Core::ElectroThermalCoupling>());
        problem->setup();
    }
};

TEST_F(Coupled2DValidationTest, CompareAgainstVtuResult) {
    // FIX: Input voltage reduced to a physically realistic value.
    constexpr double V_in = 0.1; // Was 0.5
    constexpr double T_sink = 293.15;
    constexpr double bar_width = 0.02;
    constexpr double bar_height = 0.01;
    constexpr double eps = 1e-8;

    auto *emag_field = problem->getField("Voltage");
    auto *heat_field = problem->getField("Temperature");
    auto &dof_manager = problem->getDofManager();
    const auto &mesh_ref = problem->getMesh();
    ASSERT_NE(emag_field, nullptr);
    ASSERT_NE(heat_field, nullptr);

    // Track constrained nodes and edges to avoid redundant BC application
    std::set<int> constrained_volt_vertex_nodes;
    std::set<std::pair<int, int>> constrained_volt_edges;
    std::set<int> constrained_temp_vertex_nodes;
    std::set<std::pair<int, int>> constrained_temp_edges;

    // Apply boundary conditions to both vertex and edge nodes
    for (const auto &elem_ptr : mesh_ref.getElements()) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
        if (!tri_elem) continue; // Should only have TriElements in this 2D mesh

        const auto& nodes = tri_elem->getNodes();
        // Iterate through edges (0-1, 1-2, 2-0)
        for (size_t i = 0; i < nodes.size(); ++i) {
            size_t j = (i + 1) % nodes.size();
            Core::Node* node1 = nodes[i];
            Core::Node* node2 = nodes[j];

            const auto& coords1 = node1->getCoords();
            const auto& coords2 = node2->getCoords();

            // Voltage BCs
            bool edge_on_left_boundary = (std::abs(coords1[0] - 0.0) < eps && std::abs(coords2[0] - 0.0) < eps);
            bool edge_on_right_boundary = (std::abs(coords1[0] - bar_width) < eps && std::abs(coords2[0] - bar_width) < eps);

            if (edge_on_left_boundary || edge_on_right_boundary) {
                // Apply to vertex nodes
                if (constrained_volt_vertex_nodes.find(node1->getId()) == constrained_volt_vertex_nodes.end()) {
                    emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node1->getId(), "Voltage", Eigen::Vector<double, 1>(edge_on_left_boundary ? V_in : 0.0)));
                    constrained_volt_vertex_nodes.insert(node1->getId());
                }
                if (constrained_volt_vertex_nodes.find(node2->getId()) == constrained_volt_vertex_nodes.end()) {
                    emag_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node2->getId(), "Voltage", Eigen::Vector<double, 1>(edge_on_left_boundary ? V_in : 0.0)));
                    constrained_volt_vertex_nodes.insert(node2->getId());
                }

                // Apply to edge midpoint DOF if order > 1
                if (emag_field->getElementOrder() > 1) {
                    std::vector<int> edge_node_ids = {node1->getId(), node2->getId()};
                    std::sort(edge_node_ids.begin(), edge_node_ids.end());
                    std::pair<int, int> edge_key = {edge_node_ids[0], edge_node_ids[1]};

                    if (constrained_volt_edges.find(edge_key) == constrained_volt_edges.end()) {
                        int edge_dof_idx = dof_manager.getEdgeEquationIndex(edge_node_ids, "Voltage");
                        if(edge_dof_idx != -1) {
                            emag_field->addBC(std::make_unique<Core::DirichletBC>(edge_dof_idx, Eigen::Vector<double, 1>(edge_on_left_boundary ? V_in : 0.0)));
                            constrained_volt_edges.insert(edge_key);
                        }
                    }
                }
            }

            // Heat BCs - Fixed temperature on all outer surfaces (T_sink)
            bool edge_on_heat_boundary = (std::abs(coords1[0] - 0.0) < eps && std::abs(coords2[0] - 0.0) < eps) ||
                                         (std::abs(coords1[0] - bar_width) < eps && std::abs(coords2[0] - bar_width) < eps) ||
                                         (std::abs(coords1[1] - 0.0) < eps && std::abs(coords2[1] - 0.0) < eps) ||
                                         (std::abs(coords1[1] - bar_height) < eps && std::abs(coords2[1] - bar_height) < eps);

            if (edge_on_heat_boundary) {
                // Apply to vertex nodes
                if (constrained_temp_vertex_nodes.find(node1->getId()) == constrained_temp_vertex_nodes.end()) {
                    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node1->getId(), "Temperature", Eigen::Vector<double, 1>(T_sink)));
                    constrained_temp_vertex_nodes.insert(node1->getId());
                }
                if (constrained_temp_vertex_nodes.find(node2->getId()) == constrained_temp_vertex_nodes.end()) {
                    heat_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node2->getId(), "Temperature", Eigen::Vector<double, 1>(T_sink)));
                    constrained_temp_vertex_nodes.insert(node2->getId());
                }

                // Apply to edge midpoint DOF if order > 1
                if (heat_field->getElementOrder() > 1) {
                    std::vector<int> edge_node_ids = {node1->getId(), node2->getId()};
                    std::sort(edge_node_ids.begin(), edge_node_ids.end());
                    std::pair<int, int> edge_key = {edge_node_ids[0], edge_node_ids[1]};

                    if (constrained_temp_edges.find(edge_key) == constrained_temp_edges.end()) {
                        int edge_dof_idx = dof_manager.getEdgeEquationIndex(edge_node_ids, "Temperature");
                        if(edge_dof_idx != -1) {
                            heat_field->addBC(std::make_unique<Core::DirichletBC>(edge_dof_idx, Eigen::Vector<double, 1>(T_sink)));
                            constrained_temp_edges.insert(edge_key);
                        }
                    }
                }
            }
        }
    }

    // --- Solve ---
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("results.vtk");

    // --- Load Reference Data (Coordinates and Values) ---
    auto& logger = Utils::Logger::instance();
    std::vector<std::string> data_names = {"&#x6e29;&#x5ea6;", "&#x7535;&#x52bf;"};
    IO::VtuData reference_data;
    try {
        reference_data = IO::Importer::read_vtu_points_and_data(vtu_filename, data_names);
    } catch (const Exception::FileIOException& e) {
        logger.error("Could not open reference VTU file: ", vtu_filename);
        logger.error("Please ensure the COMSOL reference file is present at the correct path.");
        FAIL() << "Reference VTU file not found.";
    }


    // --- Robust Comparison via Nearest-Neighbor Search ---
    const auto& fem_temp_solution = heat_field->getSolution();
    const auto& fem_volt_solution = emag_field->getSolution();

    ASSERT_TRUE(reference_data.point_data.count("&#x6e29;&#x5ea6;"));
    ASSERT_TRUE(reference_data.point_data.count("&#x7535;&#x52bf;"));


    const auto& ref_temp_values = reference_data.point_data.at("&#x6e29;&#x5ea6;");
    const auto& ref_volt_values = reference_data.point_data.at("&#x7535;&#x52bf;");

    double max_temp_diff = 0.0;
    double max_volt_diff = 0.0;
    int matched_nodes = 0;

    // For each node in our simulation, find the closest point in the reference data
    for (const auto& sim_node : mesh_ref.getNodes()) {
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

        // Only compare if the nodes are effectively at the same location
        if (min_dist_sq < 1e-12) {
            matched_nodes++;
            int temp_dof = dof_manager.getEquationIndex(sim_node->getId(), "Temperature");
            if (temp_dof != -1) {
                max_temp_diff = std::max(max_temp_diff, std::abs(fem_temp_solution(temp_dof) - ref_temp_values[closest_ref_idx]));
            }

            int volt_dof = dof_manager.getEquationIndex(sim_node->getId(), "Voltage");
            if (volt_dof != -1) {
                max_volt_diff = std::max(max_volt_diff, std::abs(fem_volt_solution(volt_dof) - ref_volt_values[closest_ref_idx]));
            }
        }
    }

    logger.info("Successfully matched and compared ", matched_nodes, " nodes based on coordinates.");
    if(matched_nodes < mesh_ref.getNodes().size()){
        logger.warn("Could not find a match for ", mesh_ref.getNodes().size() - matched_nodes, " nodes from the simulation mesh in the reference VTU file.");
    }

    logger.info("Maximum temperature difference: ", max_temp_diff, " K");
    logger.info("Maximum voltage difference: ", max_volt_diff, " V");

    // Assert that the maximum differences are within an acceptable tolerance
    ASSERT_LT(max_temp_diff, 1e-5);
    ASSERT_LT(max_volt_diff, 1e-5);
}
