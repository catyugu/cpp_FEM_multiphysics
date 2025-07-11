#include "core/Problem.hpp"
#include "core/LinearSolver.hpp"
#include "core/Mesh.hpp"
#include "core/DOFManager.hpp"
#include "physics/PhysicsField.hpp"
#include "physics/Current1D.hpp"
#include "physics/Heat1D.hpp"
#include "physics/Current2D.hpp"
#include "physics/Heat2D.hpp"
#include "utils/SimpleLogger.hpp"
#include "io/Exporter.hpp"
#include <iostream>
#include <iomanip>

namespace Core {

// Constructor, destructor, addField, getters, and setters remain the same...
Problem::Problem(std::unique_ptr<Mesh> mesh) : mesh_(std::move(mesh)) {
    dof_manager_ = std::make_unique<DOFManager>(*mesh_);
}
Problem::~Problem() = default;
void Problem::addField(std::unique_ptr<Physics::PhysicsField> field) {
    dof_manager_->registerVariable(field->getVariableName());
    fields_.push_back(std::move(field));
}
Physics::PhysicsField* Problem::getField(const std::string& var_name) const {
    for (const auto& field : fields_) {
        if (field->getVariableName() == var_name) return field.get();
    }
    return nullptr;
}
const Mesh& Problem::getMesh() const { return *mesh_; }
const DOFManager& Problem::getDofManager() const { return *dof_manager_; }
void Problem::setup() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("--- Problem Setup ---");
    dof_manager_->build();
    for (const auto& field : fields_) {
        field->setup(*mesh_, *dof_manager_);
    }
    logger.info("--- Problem Setup Complete ---");
}
void Problem::exportResults(const std::string& filename) const {
    IO::Exporter::write_vtk(filename, *this);
}
void Problem::setIterativeSolverParameters(int max_iter, double tol) {
    max_iterations_ = max_iter;
    convergence_tolerance_ = tol;
}
void Problem::setTimeSteppingParameters(double time_step, double total_time) {
    time_step_ = time_step;
    total_time_ = total_time;
}

// --- Corrected Steady-State Solver ---
void Problem::solveSteadyState() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("\n--- Starting Steady-State Solve ---");

    auto* emag_field = dynamic_cast<Physics::Current2D*>(getField("Voltage"));
    auto* heat_field = dynamic_cast<Physics::Heat2D*>(getField("Temperature"));

    if (emag_field && heat_field) {
        logger.info("\n--- Solving Coupled 2D Electro-Thermal Problem ---");

        // 1. Ensure the heat field is unlinked before the first EMag solve.
        //    This forces EMag2D::assemble to use a default temperature, preventing a crash.
        emag_field->setCoupledHeatField(nullptr);

        // 2. Solve the EMag field first.
        logger.info("--> Solving EMag Field (assuming constant properties)...");
        emag_field->assemble();
        emag_field->applyBCs();
        LinearSolver::solve(emag_field->getStiffnessMatrix(), emag_field->getRHSVector(), emag_field->getSolution());

        // 3. Calculate Joule heat based on the solved voltage distribution
        //    and set it as the source for the thermal problem.
        auto joule_heat = emag_field->calculateJouleHeat();
        heat_field->setVolumetricHeatSource(joule_heat);

        // 4. Solve the Heat field with the calculated heat source.
        logger.info("--> Solving Heat Field...");
        heat_field->assemble();
        heat_field->applyBCs();
        LinearSolver::solve(heat_field->getStiffnessMatrix(), heat_field->getRHSVector(), heat_field->getSolution());

    } else {
        // Fallback for single-field steady-state problems
        logger.info("\n--- Solving Uncoupled Steady-State Physics ---");
        for (const auto& field : fields_) {
            logger.info("Solving for field: ", field->getName());
            field->assemble();
            field->applyBCs();
            LinearSolver::solve(field->getStiffnessMatrix(), field->getRHSVector(), field->getSolution());
        }
    }

    logger.info("\n--- Steady-State Solve Finished ---");
}
// --- Fully Updated Transient Solver ---
void Problem::solveTransient() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("\n--- Starting Transient Solve ---");
    logger.info("Time Step: ", time_step_, "s, Total Time: ", total_time_, "s");

    auto* emag_field_base = getField("Voltage");
    auto* heat_field_base = getField("Temperature");

    if (emag_field_base && heat_field_base) {
        logger.info("\n--- Solving Coupled Transient Problem ---");

        // --- FIX: Implement the full coupled transient loop ---
        if (auto* emag1d = dynamic_cast<Physics::Current1D*>(emag_field_base)) {
            emag1d->setCoupledHeatField(heat_field_base);
        } else if (auto* emag2d = dynamic_cast<Physics::Current2D*>(emag_field_base)) {
            emag2d->setCoupledHeatField(heat_field_base);
        }

        int num_steps = static_cast<int>(total_time_ / time_step_);
        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i+1)*time_step_, "s");

            emag_field_base->assemble();
            emag_field_base->applyBCs();
            LinearSolver::solve(emag_field_base->getStiffnessMatrix(), emag_field_base->getRHSVector(), emag_field_base->getSolution());

            std::vector<double> joule_heat;
            if (auto* emag1d = dynamic_cast<Physics::Current1D*>(emag_field_base)) {
                joule_heat = emag1d->calculateJouleHeat();
            } else if (auto* emag2d = dynamic_cast<Physics::Current2D*>(emag_field_base)) {
                joule_heat = emag2d->calculateJouleHeat();
            }

            if (auto* heat1d = dynamic_cast<Physics::Heat1D*>(heat_field_base)) {
                heat1d->setVolumetricHeatSource(joule_heat);
            } else if (auto* heat2d = dynamic_cast<Physics::Heat2D*>(heat_field_base)) {
                heat2d->setVolumetricHeatSource(joule_heat);
            }

            heat_field_base->assemble();
            Eigen::SparseMatrix<double> A = (heat_field_base->getMassMatrix() / time_step_) + heat_field_base->getStiffnessMatrix();
            Eigen::VectorXd b = heat_field_base->getRHSVector() + (heat_field_base->getMassMatrix() / time_step_) * heat_field_base->getPreviousSolution();
            auto A_bc = A;
            auto b_bc = b;
            for(const auto& bc : heat_field_base->getBCs()) {
                bc->apply(A_bc, b_bc);
            }
            LinearSolver::solve(A_bc, b_bc, heat_field_base->getSolution());
            heat_field_base->updatePreviousSolution();
        }

    } else if (fields_.size() == 1) {
        logger.info("\n--- Solving Single-Field Transient Problem ---");
        auto* field = fields_[0].get();
        field->assemble();

        Eigen::SparseMatrix<double> A_eff;
        Eigen::VectorXd b_eff;

        int num_steps = static_cast<int>(total_time_ / time_step_);
        for (int i = 0; i < num_steps; ++i) {
            logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i+1)*time_step_, "s");
            A_eff = (field->getMassMatrix() / time_step_) + field->getStiffnessMatrix();
            b_eff = field->getRHSVector() + (field->getMassMatrix() / time_step_) * field->getPreviousSolution();

            for (const auto& bc : field->getBCs()) {
                bc->apply(A_eff, b_eff);
            }

            LinearSolver::solve(A_eff, b_eff, field->getSolution());
            field->updatePreviousSolution();
        }
    } else {
        logger.error("Transient solver only supports single field or coupled EMag-Heat problems.");
    }

    logger.info("\n--- Transient Solve Finished ---");
    }

} // namespace Core
