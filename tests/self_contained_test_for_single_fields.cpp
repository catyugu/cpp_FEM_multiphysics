#include <gtest/gtest.h>
#include <memory>
#include <set>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "utils/SimpleLogger.hpp"
#include <cmath>

// All 1D Physics
#include "physics/Current1D.hpp"
#include "physics/Heat1D.hpp"
#include "physics/Magnetic1D.hpp"

// All 2D Physics
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "physics/Magnetic2D.hpp"

// All 3D Physics
#include "physics/Current3D.hpp"
#include "physics/Heat3D.hpp"
#include "physics/Magnetic3D.hpp"


// Common validation logic for a simple 1D problem where U(x) = x
void validate_1D_solution(Core::Problem& problem, const std::string& var_name) {
    auto* field = problem.getField(var_name);
    const auto& dof_manager = problem.getDofManager();
    const auto& solution = field->getSolution();

    // Validate vertex nodes
    for (const auto& node : problem.getMesh().getNodes()) {
        double x = node->getCoords()[0];
        int dof_idx = dof_manager.getEquationIndex(node->getId(), var_name);
        if (dof_idx != -1) {
            ASSERT_NEAR(solution(dof_idx), x, 1e-9);
        }
    }

    // Validate higher-order edge nodes
    for (const auto& elem : problem.getMesh().getElements()) {
        auto* node1 = elem->getNodes()[0];
        auto* node2 = elem->getNodes()[1];
        int edge_dof = dof_manager.getEdgeEquationIndex({node1->getId(), node2->getId()}, var_name);
        if (edge_dof != -1) {
            double x_mid = (node1->getCoords()[0] + node2->getCoords()[0]) / 2.0;
            ASSERT_NEAR(solution(edge_dof), x_mid, 1e-9);
        }
    }
}

// =================================================================
// ==================== 1D HIGHER-ORDER TESTS ======================
// =================================================================

TEST(HigherOrderSingleFieldTest, Current1D_Order2) {
    Core::Material material("Test");
    material.setProperty("electrical_conductivity", 1.0);
    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_1d_mesh(1.0, 10)));

    auto* field = new Physics::Current1D(material);
    field->setElementOrder(2);
    problem->addField(std::unique_ptr<Physics::PhysicsField>(field));
    problem->setup();

    field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), 0, "Voltage", Eigen::Vector<double, 1>(0.0)));
    field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), 10, "Voltage", Eigen::Vector<double, 1>(1.0)));

    ASSERT_NO_THROW(problem->solveSteadyState());
    validate_1D_solution(*problem, "Voltage");
}

TEST(HigherOrderSingleFieldTest, Heat1D_Order2) {
    Core::Material material("Test");
    material.setProperty("thermal_conductivity", 1.0);
    material.setProperty("specific_heat", 1.0);
    material.setProperty("density", 1.0);
    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_1d_mesh(1.0, 10)));

    auto* field = new Physics::Heat1D(material);
    field->setElementOrder(2);
    problem->addField(std::unique_ptr<Physics::PhysicsField>(field));
    problem->setup();

    field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), 0, "Temperature", Eigen::Vector<double, 1>(0.0)));
    field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), 10, "Temperature", Eigen::Vector<double, 1>(1.0)));

    ASSERT_NO_THROW(problem->solveSteadyState());
    validate_1D_solution(*problem, "Temperature");
}

TEST(HigherOrderSingleFieldTest, Magnetic1D_Order2) {
    Core::Material material("Test");
    material.setProperty("magnetic_permeability", 1.0);
    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_1d_mesh(1.0, 10)));

    auto* field = new Physics::Magnetic1D(material);
    field->setElementOrder(2);
    problem->addField(std::unique_ptr<Physics::PhysicsField>(field));
    problem->setup();

    field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), 0, "MagneticPotential", Eigen::Vector<double, 1>(0.0)));
    field->addBC(std::make_unique<Core::DirichletBC>(problem->getDofManager(), 10, "MagneticPotential", Eigen::Vector<double, 1>(1.0)));

    ASSERT_NO_THROW(problem->solveSteadyState());
    validate_1D_solution(*problem, "MagneticPotential");
}

// =================================================================
// ==================== 2D HIGHER-ORDER TESTS ======================
// =================================================================

// Common setup and validation for a 2D problem where U(x,y) = x
void setup_and_validate_2D_problem(const std::string& field_name, Physics::PhysicsField* field) {
    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_2d_mesh(1.0, 1.0, 5, 5)));
    field->setElementOrder(2);
    problem->addField(std::unique_ptr<Physics::PhysicsField>(field));
    problem->setup();

    const auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();
    std::set<int> constrained_vertex_nodes;
    std::set<std::pair<int, int>> constrained_edges;

    for (const auto& elem : mesh_ref.getElements()) {
        const auto& nodes = elem->getNodes();
        for (size_t i = 0; i < nodes.size(); ++i) {
            size_t j = (i + 1) % nodes.size();
            Core::Node* node1 = nodes[i];
            Core::Node* node2 = nodes[j];
            double x1 = node1->getCoords()[0], x2 = node2->getCoords()[0];

            if ((std::abs(x1 - 0.0) < 1e-9 && std::abs(x2 - 0.0) < 1e-9) || (std::abs(x1 - 1.0) < 1e-9 && std::abs(x2 - 1.0) < 1e-9)) {
                if (constrained_vertex_nodes.find(node1->getId()) == constrained_vertex_nodes.end()) {
                    field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node1->getId(), field_name, Eigen::Vector<double, 1>(x1)));
                    constrained_vertex_nodes.insert(node1->getId());
                }
                if (constrained_vertex_nodes.find(node2->getId()) == constrained_vertex_nodes.end()) {
                    field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node2->getId(), field_name, Eigen::Vector<double, 1>(x2)));
                    constrained_vertex_nodes.insert(node2->getId());
                }
                std::vector<int> edge_ids = {node1->getId(), node2->getId()};
                std::sort(edge_ids.begin(), edge_ids.end());
                std::pair<int, int> edge_key = {edge_ids[0], edge_ids[1]};
                if (constrained_edges.find(edge_key) == constrained_edges.end()) {
                    int edge_dof = dof_manager.getEdgeEquationIndex(edge_ids, field_name);
                    if(edge_dof != -1) field->addBC(std::make_unique<Core::DirichletBC>(edge_dof, Eigen::Vector<double, 1>((x1 + x2) / 2.0)));
                    constrained_edges.insert(edge_key);
                }
            }
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());

    const auto& solution = field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        ASSERT_NEAR(solution(dof_manager.getEquationIndex(node->getId(), field_name)), node->getCoords()[0], 1e-9);
    }
}


TEST(HigherOrderSingleFieldTest, Current2D_Order2) {
    Core::Material material("Test");
    material.setProperty("electrical_conductivity", 1.0);
    setup_and_validate_2D_problem("Voltage", new Physics::Current2D(material));
}

TEST(HigherOrderSingleFieldTest, Heat2D_Order2) {
    Core::Material material("Test");
    material.setProperty("thermal_conductivity", 1.0);
    material.setProperty("density", 1.0);
    material.setProperty("specific_heat", 1.0);
    setup_and_validate_2D_problem("Temperature", new Physics::Heat2D(material));
}

TEST(HigherOrderSingleFieldTest, Magnetic2D_Order2) {
    Core::Material material("Test");
    material.setProperty("magnetic_permeability", 1.0);
    setup_and_validate_2D_problem("MagneticPotential", new Physics::Magnetic2D(material));
}


// =================================================================
// ==================== 3D HIGHER-ORDER TESTS ======================
// =================================================================
// NOTE: Re-using the test case from the previous issue, which is now correct.
TEST(HigherOrderSingleFieldTest, Magnetic3D_Order2) {
    Core::Material material("TestMaterial");
    material.setProperty("magnetic_permeability", 1.0);
    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(1.0, 1.0, 1.0, 2, 2, 2)));
    auto* magnetic_field = new Physics::Magnetic3D(material);
    magnetic_field->setElementOrder(2);
    problem->addField(std::unique_ptr<Physics::PhysicsField>(magnetic_field));
    problem->setup();
    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();

    std::set<int> constrained_vertex_nodes;
    std::set<std::pair<int, int>> constrained_edges;

    for(const auto& elem : mesh_ref.getElements()) {
        const auto& nodes = elem->getNodes();
        const std::vector<std::vector<int>> faces = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}};
        for(const auto& face : faces) {
            Core::Node* n1 = nodes[face[0]], *n2 = nodes[face[1]], *n3 = nodes[face[2]];
            double x1 = n1->getCoords()[0], x2 = n2->getCoords()[0], x3 = n3->getCoords()[0];

            if((std::abs(x1-0.0)<1e-9 && std::abs(x2-0.0)<1e-9 && std::abs(x3-0.0)<1e-9) ||
               (std::abs(x1-1.0)<1e-9 && std::abs(x2-1.0)<1e-9 && std::abs(x3-1.0)<1e-9)) {

                const std::vector<std::pair<int,int>> face_edges = {{0,1}, {1,2}, {2,0}};
                for(const auto& edge_indices : face_edges) {
                    Core::Node* node_a = nodes[face[edge_indices.first]];
                    Core::Node* node_b = nodes[face[edge_indices.second]];

                    if (constrained_vertex_nodes.find(node_a->getId()) == constrained_vertex_nodes.end()) {
                        magnetic_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node_a->getId(), "MagneticPotential", Eigen::Vector<double, 1>(node_a->getCoords()[0])));
                        constrained_vertex_nodes.insert(node_a->getId());
                    }
                    if (constrained_vertex_nodes.find(node_b->getId()) == constrained_vertex_nodes.end()) {
                        magnetic_field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node_b->getId(), "MagneticPotential", Eigen::Vector<double, 1>(node_b->getCoords()[0])));
                        constrained_vertex_nodes.insert(node_b->getId());
                    }

                    std::vector<int> edge_dof_key = {node_a->getId(), node_b->getId()};
                    std::sort(edge_dof_key.begin(), edge_dof_key.end());
                    std::pair<int, int> edge_key_pair = {edge_dof_key[0], edge_dof_key[1]};
                    if (constrained_edges.find(edge_key_pair) == constrained_edges.end()) {
                        int edge_dof_idx = dof_manager.getEdgeEquationIndex(edge_dof_key, "MagneticPotential");
                        if (edge_dof_idx != -1) {
                            magnetic_field->addBC(std::make_unique<Core::DirichletBC>(edge_dof_idx, Eigen::Vector<double, 1>((node_a->getCoords()[0] + node_b->getCoords()[0]) / 2.0)));
                        }
                        constrained_edges.insert(edge_key_pair);
                    }
                }
            }
        }
    }
    ASSERT_NO_THROW(problem->solveSteadyState());
}