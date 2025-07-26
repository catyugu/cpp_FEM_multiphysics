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
    material.setProperty("thermal_capacity", 1.0);
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
    material.setProperty("thermal_capacity", 1.0);
    setup_and_validate_2D_problem("Temperature", new Physics::Heat2D(material));
}


// =================================================================
// ==================== 3D HIGHER-ORDER TESTS ======================
// =================================================================

// ** NEW ** Common setup and validation for a 3D problem where U(x,y,z) = x
// Common setup and validation for a 3D problem where U(x,y,z) = x
void setup_and_validate_3D_problem(const std::string& field_name, Physics::PhysicsField* field) {
    auto problem = std::make_unique<Core::Problem>(std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(1.0, 1.0, 1.0, 5, 5, 5)));

    field->setElementOrder(2); // Assuming Order 2 for this test
    problem->addField(std::unique_ptr<Physics::PhysicsField>(field));
    problem->setup();

    auto& dof_manager = problem->getDofManager();
    const auto& mesh_ref = problem->getMesh();
    constexpr double eps = 1e-9;

    // Apply boundary conditions to all NODES on the exterior faces. This part was correct.
    for (const auto& node : mesh_ref.getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < eps || std::abs(coords[0] - 1.0) < eps ||
            std::abs(coords[1] - 0.0) < eps || std::abs(coords[1] - 1.0) < eps ||
            std::abs(coords[2] - 0.0) < eps || std::abs(coords[2] - 1.0) < eps) {
            field->addBC(std::make_unique<Core::DirichletBC>(dof_manager, node->getId(), field_name, Eigen::Vector<double, 1>(coords[0])));
        }
    }

    // Apply boundary conditions to all higher-order EDGE DOFs on the exterior faces.
    for (const auto& elem : mesh_ref.getElements()) {
        auto element_nodes = elem->getNodes();
        for (size_t i = 0; i < element_nodes.size(); ++i) {
            for (size_t j = i + 1; j < element_nodes.size(); ++j) {
                 auto* node1 = element_nodes[i];
                 auto* node2 = element_nodes[j];
                 const auto& coords1 = node1->getCoords();
                 const auto& coords2 = node2->getCoords();

                 bool node1_on_boundary = (std::abs(coords1[0] - 0.0) < eps || std::abs(coords1[0] - 1.0) < eps || std::abs(coords1[1] - 0.0) < eps || std::abs(coords1[1] - 1.0) < eps || std::abs(coords1[2] - 0.0) < eps || std::abs(coords1[2] - 1.0) < eps);
                 bool node2_on_boundary = (std::abs(coords2[0] - 0.0) < eps || std::abs(coords2[0] - 1.0) < eps || std::abs(coords2[1] - 0.0) < eps || std::abs(coords2[1] - 1.0) < eps || std::abs(coords2[2] - 0.0) < eps || std::abs(coords2[2] - 1.0) < eps);

                 if (node1_on_boundary && node2_on_boundary) {
                    // --- FIX START ---
                    // This new check ensures the edge itself lies flat on a boundary face, not just its endpoints.
                    // It does this by checking if the nodes share a common x, y, or z boundary coordinate.
                    bool edge_is_truly_on_boundary =
                        (std::abs(coords1[0] - coords2[0]) < eps && (std::abs(coords1[0] - 0.0) < eps || std::abs(coords1[0] - 1.0) < eps)) ||
                        (std::abs(coords1[1] - coords2[1]) < eps && (std::abs(coords1[1] - 0.0) < eps || std::abs(coords1[1] - 1.0) < eps)) ||
                        (std::abs(coords1[2] - coords2[2]) < eps && (std::abs(coords1[2] - 0.0) < eps || std::abs(coords1[2] - 1.0) < eps));

                    if (edge_is_truly_on_boundary) {
                        std::vector<int> edge_dof_key = {node1->getId(), node2->getId()};
                        std::sort(edge_dof_key.begin(), edge_dof_key.end());
                        int edge_dof_idx = dof_manager.getEdgeEquationIndex(edge_dof_key, field_name);

                        if (edge_dof_idx != -1) {
                            double x_mid = (coords1[0] + coords2[0]) / 2.0;
                            field->addBC(std::make_unique<Core::DirichletBC>(edge_dof_idx, Eigen::Vector<double, 1>(x_mid)));
                        }
                    }
                    // --- FIX END ---
                 }
            }
        }
    }

    ASSERT_NO_THROW(problem->solveSteadyState());

    // Validate the solution U(x,y,z) = x
    const auto& solution = field->getSolution();
    for (const auto& node : mesh_ref.getNodes()) {
        int dof_idx = dof_manager.getEquationIndex(node->getId(), field_name);
        if (dof_idx != -1) {
            ASSERT_NEAR(solution(dof_idx), node->getCoords()[0], 1e-9);
        }
    }
     for(const auto& elem : mesh_ref.getElements()){
        auto element_nodes = elem->getNodes();
        for(size_t i = 0; i < element_nodes.size(); ++i){
            for(size_t j = i + 1; j < element_nodes.size(); ++j){
                 std::vector<int> edge_dof_key = {element_nodes[i]->getId(), element_nodes[j]->getId()};
                 std::sort(edge_dof_key.begin(), edge_dof_key.end());
                 int edge_dof_idx = dof_manager.getEdgeEquationIndex(edge_dof_key, field_name);
                 if(edge_dof_idx != -1){
                     double x_mid = (element_nodes[i]->getCoords()[0] + element_nodes[j]->getCoords()[0]) / 2.0;
                     ASSERT_NEAR(solution(edge_dof_idx), x_mid, 1e-9);
                 }
            }
        }
     }
}

// ** NEW ** Test for 3D Heat Conduction
TEST(HigherOrderSingleFieldTest, Heat3D_Order2) {
    Core::Material material("TestMaterial");
    material.setProperty("thermal_conductivity", 1.0);
    material.setProperty("density", 1.0);
    material.setProperty("thermal_capacity", 1.0);
    setup_and_validate_3D_problem("Temperature", new Physics::Heat3D(material));
}

// ** NEW ** Test for 3D Current Conduction
TEST(HigherOrderSingleFieldTest, Current3D_Order2) {
    Core::Material material("TestMaterial");
    material.setProperty("electrical_conductivity", 1.0);
    setup_and_validate_3D_problem("Voltage", new Physics::Current3D(material));
}

