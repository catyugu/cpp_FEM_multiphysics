#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "mesh.h"
#include "Solver.h"
#include "SimpleLogger.hpp"
#include <Eigen/Sparse>
#include <vector>
#include <memory>
#include <functional>

extern SimpleLogger::Logger& logger;

// A temporary struct for params, can be expanded for other physics
struct Params {
    double k{397.0};   // Thermal conductivity of copper
    double rou{8940.0}; // Density of copper
    double Cp{385.0};  // Specific heat capacity of copper
};

// The main analysis class, replacing the old HeatField class.
class Analysis {
public:
    explicit Analysis(std::unique_ptr<Mesh> mesh_ptr);
    
    // Public API for setting up and running the simulation
    void initialize_physics();
    void assemble_initial_condition(const std::function<double(const Eigen::Vector3d& coord)>& init_func);
    void assemble_boundary_condition(const std::function<double(const Eigen::Vector3d& coord)>& border_func);
    void run_simulation(int num_steps, double time_step);
    void write_vtk(const std::string& filename) const;

private:
    std::unique_ptr<Mesh> mesh;
    std::vector<std::unique_ptr<FiniteElement>> elements;
    Params params; // For this example, params are stored here

    // Global system matrices and vectors
    Eigen::SparseMatrix<double> global_stiffness;
    Eigen::SparseMatrix<double> global_mass;
    Eigen::VectorXd solution;
    Eigen::VectorXd boundary_values;
    std::vector<bool> is_border_node;

    // Private methods
    void assemble_global_matrices();
    void identify_border_nodes();
    void solve_step(double time_step);
};

#endif //ANALYSIS_H