#ifndef HEAT_FIELD_H
#define HEAT_FIELD_H

#include "SimpleLogger.hpp"
#include "mesh.h"
#include "functional"
#include <memory>
#include <string>
#include <Eigen/SparseLU> // <--- Include the Sparse header

extern SimpleLogger::Logger& logger;

struct HeatParams {
    double k{300.0};   // Thermal conductivity of copper
    double rou{8940.0}; // Density of copper
    double Cp{385.0};  // Specific heat capacity of copper
};

class HeatField {
public:
    std::unique_ptr<Mesh> mesh;
    Eigen::VectorXd init_condition;
    Eigen::VectorXd border_condition;

    // --- CHANGE TO SPARSE MATRICES ---
    Eigen::SparseMatrix<double> stiffness_matrix;
    Eigen::SparseMatrix<double> mass_matrix;
    // ---------------------------------

    Eigen::VectorXd F;
    Eigen::VectorXd heat;

    std::function<double(const Eigen::Vector3d& coord)> force_heat;
    HeatParams params;

    explicit HeatField(std::unique_ptr<Mesh> mesh_ptr);

    void apply_init_condition(const std::function<double(const Eigen::Vector3d& coord)>& init_func);
    void apply_border_condition(const std::function<double(const Eigen::Vector3d& coord)>& border_func);
    void apply_extern_field(const std::function<double(const Eigen::Vector3d& coord)>& source_func);
    void step_forward(double time_step);
    void write_vtk(const std::string& filename) const;

private:
    std::vector<bool> is_border_node;
    bool are_matrices_assembled{false};

    void identify_border_nodes();
    void assemble_matrices();
};

#endif //HEAT_FIELD_H