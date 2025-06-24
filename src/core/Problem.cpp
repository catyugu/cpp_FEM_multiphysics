#include "core/Problem.hpp"
#include "core/LinearSolver.hpp"
#include "core/Mesh.hpp"
#include "core/DOFManager.hpp"
#include "physics/PhysicsField.hpp"
#include "physics/EMag1D.hpp"
#include "physics/Heat1D.hpp"
#include "physics/EMag2D.hpp"
#include "physics/Heat2D.hpp"
#include "utils/SimpleLogger.hpp"
#include "io/Exporter.hpp"
#include <iostream>
#include <iomanip>

namespace Core {

// Constructor, destructor, addField, getters... (remain unchanged)
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

void Problem::solveSteadyState() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("\n--- Starting Steady-State Solve ---");

    auto* emag2d = dynamic_cast<Physics::EMag2D*>(getField("Voltage"));
    auto* heat2d = dynamic_cast<Physics::Heat2D*>(getField("Temperature"));

    if (emag2d && heat2d) {
        // --- Strongly Coupled 2D Iterative Solve ---
        logger.info("\n--- Solving Strongly Coupled 2D Electro-Thermal Problem ---");
        logger.info("Max Iterations: ", max_iterations_, ", Tolerance: ", convergence_tolerance_);

        emag2d->setCoupledHeatField(heat2d);

        double T_max_prev = 0.0;
        for (int i = 0; i < max_iterations_; ++i) {
            logger.info("\n--- Iteration ", i + 1, " ---");

            emag2d->assemble();
            emag2d->applyBCs();
            LinearSolver::solve(emag2d->getStiffnessMatrix(), emag2d->getRHSVector(), emag2d->getSolution());

            auto joule_heat = emag2d->calculateJouleHeat();
            heat2d->setVolumetricHeatSource(joule_heat);

            heat2d->assemble();
            heat2d->applyBCs();
            LinearSolver::solve(heat2d->getStiffnessMatrix(), heat2d->getRHSVector(), heat2d->getSolution());

            double T_max_curr = heat2d->getSolution().maxCoeff();
            double relative_error = std::abs(T_max_curr - T_max_prev) / (T_max_curr + 1e-9);
            logger.info("Current max temperature: ", T_max_curr, " K. Relative change: ", relative_error);

            if (relative_error < convergence_tolerance_) {
                logger.info("Solution converged after ", i + 1, " iterations.");
                break;
            }
            if(i == max_iterations_ - 1) {
                logger.warn("Solver did not converge within the maximum number of iterations.");
            }
            T_max_prev = T_max_curr;
        }

    } else {
        // --- Fallback to simple uncoupled solve for all other steady-state cases ---
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

void Problem::solveTransient() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("\n--- Starting Transient Solve ---");
    logger.info("Time Step: ", time_step_, "s, Total Time: ", total_time_, "s");

    if (fields_.size() != 1) {
        logger.error("Transient solver currently only supports single-field problems.");
        return;
    }
    auto* field = fields_[0].get();

    field->assemble();

    int num_steps = static_cast<int>(total_time_ / time_step_);
    for (int i = 0; i < num_steps; ++i) {
        logger.info("Time Step ", i + 1, " / ", num_steps, ", Time = ", (i+1)*time_step_, "s");

        Eigen::SparseMatrix<double> A = (field->getMassMatrix() / time_step_) + field->getStiffnessMatrix();
        Eigen::VectorXd b = field->getRHSVector() + (field->getMassMatrix() / time_step_) * field->getPreviousSolution();

        auto A_bc = A;
        auto b_bc = b;
        for(const auto& bc : field->getBCs()) {
            bc->apply(A_bc, b_bc);
        }

        LinearSolver::solve(A_bc, b_bc, field->getSolution());

        field->updatePreviousSolution();
    }

    logger.info("\n--- Transient Solve Finished ---");
}

} // namespace Core
