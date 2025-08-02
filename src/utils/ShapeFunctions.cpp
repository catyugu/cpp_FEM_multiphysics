#include "utils/ShapeFunctions.hpp"
#include <stdexcept>
#include <vector>

namespace Utils {

// --- 1D Line Shape Functions ---
Eigen::VectorXd ShapeFunctions::getLineShapeFunctions(int order, double xi) {
    if (order == 1) { // Linear, nodes at xi = -1, 1
        Eigen::Vector2d N;
        N << (1.0 - xi) / 2.0, (1.0 + xi) / 2.0;
        return N;
    }
    if (order == 2) { // Quadratic, standard node order: xi = -1, 1, 0
        Eigen::VectorXd N(3);
        N(0) = 0.5 * xi * (xi - 1.0);  // N for node at xi = -1
        N(1) = 0.5 * xi * (xi + 1.0);  // N for node at xi = +1
        N(2) = (1.0 - xi * xi);        // N for node at xi = 0
        return N;
    }
    throw std::invalid_argument("Line shape function order > 2 not implemented with standard ordering.");
}

Eigen::VectorXd ShapeFunctions::getLineShapeFunctionDerivatives(int order, double xi) {
    if (order == 1) {
        Eigen::Vector2d dN_dxi;
        dN_dxi << -0.5, 0.5;
        return dN_dxi;
    }
    if (order == 2) {
        Eigen::VectorXd dN_dxi(3);
        dN_dxi(0) = 0.5 * (2.0 * xi - 1.0); // dN0/dxi
        dN_dxi(1) = 0.5 * (2.0 * xi + 1.0); // dN1/dxi
        dN_dxi(2) = -2.0 * xi;              // dN2/dxi
        return dN_dxi;
    }
    throw std::invalid_argument("Line shape function derivative order > 2 not implemented with standard ordering.");
}

// --- 2D Triangle Shape Functions (unchanged, but provided for completeness) ---
Eigen::VectorXd ShapeFunctions::getTriShapeFunctions(int order, double xi, double eta) {
    double L1 = 1.0 - xi - eta;
    double L2 = xi;
    double L3 = eta;

    if (order == 1) {
        Eigen::Vector3d N;
        N << L1, L2, L3;
        return N;
    }
    if (order == 2) {
        Eigen::VectorXd N(6);
        N << L1 * (2 * L1 - 1),
             L2 * (2 * L2 - 1),
             L3 * (2 * L3 - 1),
             4 * L1 * L2,
             4 * L2 * L3,
             4 * L3 * L1;
        return N;
    }
    throw std::invalid_argument("Triangle shape function order > 2 not yet implemented.");
}

Eigen::MatrixXd ShapeFunctions::getTriShapeFunctionDerivatives(int order, double xi, double eta) {
    double L1 = 1.0 - xi - eta;
    double L2 = xi;
    double L3 = eta;

    if (order == 1) {
        Eigen::MatrixXd dN(3, 2);
        dN << -1.0, -1.0,
               1.0,  0.0,
               0.0,  1.0;
        return dN;
    }
    if (order == 2) {
        Eigen::MatrixXd dN(6, 2);
        dN(0, 0) = -1.0 * (4 * L1 - 1);
        dN(1, 0) = 4 * L2 - 1;
        dN(2, 0) = 0;
        dN(3, 0) = 4 * (L1 - L2);
        dN(4, 0) = 4 * L3;
        dN(5, 0) = -4 * L3;
        dN(0, 1) = -1.0 * (4 * L1 - 1);
        dN(1, 1) = 0;
        dN(2, 1) = 4 * L3 - 1;
        dN(3, 1) = -4 * L2;
        dN(4, 1) = 4 * L2;
        dN(5, 1) = 4 * (L1 - L3);
        return dN;
    }
    throw std::invalid_argument("Triangle shape function derivative order > 2 not yet implemented.");
}


// --- 3D Tetrahedron Shape Functions (unchanged, but provided for completeness) ---
Eigen::VectorXd ShapeFunctions::getTetShapeFunctions(int order, double xi, double eta, double zeta) {
    double L1 = 1.0 - xi - eta - zeta;
    double L2 = xi;
    double L3 = eta;
    double L4 = zeta;

    if (order == 1) {
        Eigen::Vector4d N;
        N << L1, L2, L3, L4;
        return N;
    }
    if (order == 2) {
        Eigen::VectorXd N(10);
        N << L1 * (2*L1-1), L2 * (2*L2-1), L3 * (2*L3-1), L4 * (2*L4-1), // Vertices
             4*L1*L2, 4*L1*L3, 4*L1*L4, 4*L2*L3, 4*L2*L4, 4*L3*L4;      // Edges
        return N;
    }
    throw std::invalid_argument("Tetrahedron shape function order > 2 not yet implemented.");
}

Eigen::MatrixXd ShapeFunctions::getTetShapeFunctionDerivatives(int order, double xi, double eta, double zeta) {
    if (order == 1) {
        Eigen::MatrixXd dN(4, 3);
        dN.row(0) << -1, -1, -1;
        dN.row(1) <<  1,  0,  0;
        dN.row(2) <<  0,  1,  0;
        dN.row(3) <<  0,  0,  1;
        return dN;
    }
     if (order == 2) {
        double L1 = 1.0 - xi - eta - zeta, L2 = xi, L3 = eta, L4 = zeta;
        Eigen::MatrixXd dN(10, 3);
        dN.row(0) << -1*(4*L1-1), -1*(4*L1-1), -1*(4*L1-1);
        dN.row(1) <<  1*(4*L2-1),  0,           0;
        dN.row(2) <<  0,           1*(4*L3-1),  0;
        dN.row(3) <<  0,           0,           1*(4*L4-1);
        dN.row(4) <<  4*(L1-L2),  -4*L2, -4*L2;
        dN.row(5) <<  -4*L3, 4*(L1-L3), -4*L3;
        dN.row(6) <<  -4*L4, -4*L4, 4*(L1-L4);
        dN.row(7) <<   4*L3,  4*L2,  0;
        dN.row(8) <<   4*L4,   0,   4*L2;
        dN.row(9) <<   0,    4*L4,  4*L3;
        return dN;
    }
    throw std::invalid_argument("Tetrahedron shape function derivative order > 2 not yet implemented.");
}


} // namespace Utils