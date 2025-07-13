#include "core/Problem.hpp"
#include <solver/LinearSolver.hpp> // Make sure this is included for SolverType
#include <core/mesh/Mesh.hpp>
#include "core/DOFManager.hpp"
#include "physics/PhysicsField.hpp"
#include "solver/SolverFactory.hpp"
#include "io/Exporter.hpp"
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"
#include <iostream>
#include <iomanip>

namespace Core {

    Problem::Problem(std::unique_ptr<Mesh> mesh) : mesh_(std::move(mesh)) {
        dof_manager_ = std::make_unique<DOFManager>(*mesh_);
    }

    Problem::~Problem() = default;

    void Problem::addField(std::unique_ptr<Physics::PhysicsField> field) {
        coupling_manager_.registerField(*field);
        dof_manager_->registerVariable(field->getVariableName());
        fields_.push_back(std::move(field));
    }

    void Problem::setup() {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("--- Problem Setup ---");
        dof_manager_->build();
        for (const auto &field: fields_) {
            field->setup(*mesh_, *dof_manager_);
        }
        coupling_manager_.setupCouplings();
        solver_ = Solver::SolverFactory::createSolver(*this);
        logger.info("--- Problem Setup Complete ---");
    }

    void Problem::solveSteadyState() {
        if (solver_) {
            solver_->solveSteadyState(*this);
        } else {
            throw Exception::ConfigurationException("Solver not initialized. Did you forget to call setup()?");
        }
    }

    void Problem::solveTransient() {
        if (solver_) {
            solver_->solveTransient(*this);
        } else {
            throw Exception::ConfigurationException("Solver not initialized. Did you forget to call setup()?");
        }
    }

    void Problem::exportResults(const std::string &filename) const {
        IO::Exporter::write_vtk(filename, *this);
    }

    const Mesh& Problem::getMesh() const {
        return *mesh_;
    }

    const DOFManager& Problem::getDofManager() const {
        return *dof_manager_;
    }

    Physics::PhysicsField* Problem::getField(const std::string& var_name) const {
        for (const auto& field : fields_) {
            if (field->getVariableName() == var_name) {
                return field.get();
            }
        }
        return nullptr;
    }

    const std::vector<std::unique_ptr<Physics::PhysicsField>>& Problem::getFields() const {
        return fields_;
    }

    void Problem::setIterativeSolverParameters(int max_iter, double tol) {
        if (max_iter <= 0) throw std::invalid_argument("Maximum number of iterations must be positive.");
        if (tol <= 0) throw std::invalid_argument("Convergence tolerance must be positive.");
        max_iterations_ = max_iter;
        convergence_tolerance_ = tol;
    }

    void Problem::setTimeSteppingParameters(double time_step, double total_time) {
        if (time_step <= 0) throw std::invalid_argument("Time step must be positive.");
        if (total_time <= 0) throw std::invalid_argument("Total time must be positive.");
        time_step_ = time_step;
        total_time_ = total_time;
    }

    // New method implementation
    void Problem::setLinearSolverType(Solver::SolverType type) {
        linear_solver_type_ = type;
        SimpleLogger::Logger::instance().info("Problem: Set linear solver type to ",
                                             (type == Solver::SolverType::LU ? "LU" : "BiCGSTAB"));
    }

} // namespace Core