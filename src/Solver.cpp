#include "Solver.h"
#include "HeatField.h"

void HeatTetrahedron::calculate_local_matrices(const Node* const* element_nodes,
                                               const HeatParams& params,
                                               Eigen::MatrixXd& local_stiffness,
                                               Eigen::MatrixXd& local_mass) const {
    // 1. Get nodal coordinates
    const Eigen::Vector3d& p1 = element_nodes[0]->coordinates;
    const Eigen::Vector3d& p2 = element_nodes[1]->coordinates;
    const Eigen::Vector3d& p3 = element_nodes[2]->coordinates;
    const Eigen::Vector3d& p4 = element_nodes[3]->coordinates;

    // 2. Calculate the volume of the tetrahedron
    double volume = std::abs((p1 - p4).dot((p2 - p4).cross(p3 - p4))) / 6.0;
    if (volume < 1e-12) {
        local_stiffness.setZero(4, 4);
        local_mass.setZero(4, 4);
        return; // Skip degenerate elements
    }

    // 3. Calculate shape function gradients (B matrix) for 3D
    Eigen::Matrix4d M_coords;
    M_coords << 1, p1.x(), p1.y(), p1.z(),
                1, p2.x(), p2.y(), p2.z(),
                1, p3.x(), p3.y(), p3.z(),
                1, p4.x(), p4.y(), p4.z();

    Eigen::Matrix<double, 3, 4> B = M_coords.inverse().bottomRows<3>();

    // 4. Calculate local stiffness matrix: k_e = volume * k * B^T * B
    local_stiffness = params.k * volume * (B.transpose() * B);

    // 5. Calculate local mass matrix (consistent formulation)
    local_mass.setConstant(4, 4, 1.0); // Set off-diagonals to 1
    local_mass.diagonal().setConstant(2.0); // Set diagonals to 2
    local_mass *= (params.rou * params.Cp * volume / 20.0);
}