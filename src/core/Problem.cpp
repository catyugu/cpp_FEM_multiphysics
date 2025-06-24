#include "core/Problem.hpp"
#include "core/LinearSolver.hpp"
#include "core/Mesh.hpp"
#include "core/DOFManager.hpp"
#include "physics/PhysicsField.hpp"
#include "physics/EMag1D.hpp"
#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <iostream>
#include <iomanip>

namespace Core {

Problem::Problem(std::unique_ptr<Mesh> mesh) : mesh_(std::move(mesh)) {
    dof_manager_ = std::make_unique<DOFManager>(*mesh_);
}

Problem::~Problem() = default;

void Problem::addField(std::unique_ptr<Physics::PhysicsField> field) {
    dof_manager_->registerVariable(field->getVariableName());
    fields_.push_back(std::move(field));
}

Physics::PhysicsField* Problem::getField(const std::string& var_name) {
    for (const auto& field : fields_) {
        if (field->getVariableName() == var_name) {
            return field.get();
        }
    }
    return nullptr;
}

Mesh& Problem::getMesh() const { return *mesh_; }
DOFManager& Problem::getDofManager() const { return *dof_manager_; }


void Problem::setup() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("--- Problem Setup ---");

    dof_manager_->build();

    for (const auto& field : fields_) {
        field->setup(*mesh_, *dof_manager_);
    }
    logger.info("--- Problem Setup Complete ---");
}

void Problem::solve() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("\n--- Starting Problem Solve Sequence ---");

    // --- FIX: More generic solve loop ---

    // 1. Handle Coupling: Check for the specific electro-thermal case first.
    Physics::EMag1D* emag_field = dynamic_cast<Physics::EMag1D*>(getField("Voltage"));
    Physics::Heat1D* heat_field = dynamic_cast<Physics::Heat1D*>(getField("Temperature"));

    if (emag_field && heat_field) {
        // --- Coupled Solve ---
        logger.info("\n--- Solving Coupled Electro-Thermal Problem ---");
        // a. Solve EMag
        emag_field->assemble();
        emag_field->applyBCs();
        LinearSolver::solve(emag_field->getStiffnessMatrix(), emag_field->getRHSVector(), emag_field->getSolution());

        // b. Perform coupling and solve Heat
        auto joule_heat = emag_field->calculateJouleHeat();
        heat_field->setVolumetricHeatSource(joule_heat);
        heat_field->assemble();
        heat_field->applyBCs();
        LinearSolver::solve(heat_field->getStiffnessMatrix(), heat_field->getRHSVector(), heat_field->getSolution());

    } else {
        // --- Uncoupled Solve ---
        // If it's not the specific coupled case, solve each field present, one by one.
        logger.info("\n--- Solving Uncoupled Physics ---");
        for (const auto& field : fields_) {
            logger.info("Solving for field: ", field->getName());
            field->assemble();
            field->applyBCs();
            LinearSolver::solve(field->getStiffnessMatrix(), field->getRHSVector(), field->getSolution());
        }
    }

    logger.info("\n--- Problem Solve Sequence Finished ---");
}

} // namespace Core
