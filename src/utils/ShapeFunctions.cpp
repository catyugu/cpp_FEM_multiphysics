#include "utils/ShapeFunctions.hpp"
#include <stdexcept>
#include <vector>

namespace Utils {

// Helper for Lagrangian polynomials
double lagrange_basis(int order, int i, double xi, const std::vector<double>& points) {
    double result = 1.0;
    for (int j = 0; j <= order; ++j) {
        if (i == j) continue;
        result *= (xi - points[j]) / (points[i] - points[j]);
    }
    return result;
}

double lagrange_basis_derivative(int order, int i, double xi, const std::vector<double>& points) {
    double result = 0.0;
    for (int k = 0; k <= order; ++k) {
        if (k == i) continue;
        double term = 1.0 / (points[i] - points[k]);
        for (int j = 0; j <= order; ++j) {
            if (j == i || j == k) continue;
            term *= (xi - points[j]) / (points[i] - points[j]);
        }
        result += term;
    }
    return result;
}


// --- 1D Line Shape Functions ---
Eigen::VectorXd ShapeFunctions::getLineShapeFunctions(int order, double xi) {
    int num_nodes = order + 1;
    Eigen::VectorXd N(num_nodes);
    std::vector<double> points(num_nodes);
    for(int i = 0; i < num_nodes; ++i) {
        points[i] = -1.0 + 2.0 * i / order;
    }

    for (int i = 0; i < num_nodes; ++i) {
        N(i) = lagrange_basis(order, i, xi, points);
    }
    return N;
}

Eigen::VectorXd ShapeFunctions::getLineShapeFunctionDerivatives(int order, double xi) {
    int num_nodes = order + 1;
    Eigen::VectorXd dN_dxi(num_nodes);
    std::vector<double> points(num_nodes);
    for(int i = 0; i < num_nodes; ++i) {
        points[i] = -1.0 + 2.0 * i / order;
    }

    for (int i = 0; i < num_nodes; ++i) {
        dN_dxi(i) = lagrange_basis_derivative(order, i, xi, points);
    }
    return dN_dxi;
}


// --- 2D Triangle Shape Functions (in terms of L-coordinates) ---
Eigen::VectorXd ShapeFunctions::getTriShapeFunctions(int order, double xi, double eta) {
    double L1 = 1.0 - xi - eta;
    double L2 = xi;
    double L3 = eta;

    // This is a placeholder for a proper implementation based on L-coordinates
    // For simplicity, we only implement P1 and P2 for now.
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

    // Orders 3, 4, 5 would require a more general formula for the basis functions
    // which is significantly more complex. We will throw an error for now.
    throw std::invalid_argument("Triangle shape function order > 2 not yet implemented.");
}

Eigen::MatrixXd ShapeFunctions::getTriShapeFunctionDerivatives(int order, double xi, double eta) {
    double L1 = 1.0 - xi - eta;
    double L2 = xi;
    double L3 = eta;

    // Derivatives of L-coordinates w.r.t xi and eta
    // dL1/dxi = -1, dL1/deta = -1
    // dL2/dxi = 1,  dL2/deta = 0
    // dL3/dxi = 0,  dL3/deta = 1

    if (order == 1) {
        Eigen::MatrixXd dN(3, 2);
        // dN/dxi, dN/deta
        dN << -1.0, -1.0,  // N1
               1.0,  0.0,  // N2
               0.0,  1.0;  // N3
        return dN;
    }
    if (order == 2) {
        Eigen::MatrixXd dN(6, 2);
        // dN/dxi
        dN(0, 0) = -1.0 * (4 * L1 - 1);       // d(N1)/dxi
        dN(1, 0) = 4 * L2 - 1;                // d(N2)/dxi
        dN(2, 0) = 0;                         // d(N3)/dxi
        dN(3, 0) = 4 * (L1 - L2);             // d(N4)/dxi
        dN(4, 0) = 4 * L3;                    // d(N5)/dxi
        dN(5, 0) = -4 * L3;                   // d(N6)/dxi
        // dN/deta
        dN(0, 1) = -1.0 * (4 * L1 - 1);       // d(N1)/deta
        dN(1, 1) = 0;                         // d(N2)/deta
        dN(2, 1) = 4 * L3 - 1;                // d(N3)/deta
        dN(3, 1) = -4 * L2;                   // d(N4)/deta
        dN(4, 1) = 4 * (L2 - L3);             // d(N5)/deta
        dN(5, 1) = 4 * (L1 - L3);             // d(N6)/deta
        return dN;
    }
    throw std::invalid_argument("Triangle shape function derivative order > 2 not yet implemented.");
}


// --- 3D Tetrahedron Shape Functions ---

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
    // Derivatives of L-coords w.r.t natural coords (xi, eta, zeta)
    // dL1 = [-1, -1, -1]
    // dL2 = [ 1,  0,  0]
    // dL3 = [ 0,  1,  0]
    // dL4 = [ 0,  0,  1]

    if (order == 1) {
        Eigen::MatrixXd dN(4, 3);
        dN.row(0) << -1, -1, -1; // dN1
        dN.row(1) <<  1,  0,  0; // dN2
        dN.row(2) <<  0,  1,  0; // dN3
        dN.row(3) <<  0,  0,  1; // dN4
        return dN;
    }
     if (order == 2) {
        double L1 = 1.0 - xi - eta - zeta, L2 = xi, L3 = eta, L4 = zeta;
        Eigen::MatrixXd dN(10, 3);
        // Derivatives of vertex nodes
        dN.row(0) << -1*(4*L1-1), -1*(4*L1-1), -1*(4*L1-1); // dN1
        dN.row(1) <<  1*(4*L2-1),  0,           0;          // dN2
        dN.row(2) <<  0,           1*(4*L3-1),  0;          // dN3
        dN.row(3) <<  0,           0,           1*(4*L4-1); // dN4
        // Derivatives of edge nodes
        dN.row(4) <<  4*(L1-L2),  -4*L2, -4*L2; // dN5 (L1,L2)
        dN.row(5) <<  -4*L3, 4*(L1-L3), -4*L3; // dN6 (L1,L3)
        dN.row(6) <<  -4*L4, -4*L4, 4*(L1-L4); // dN7 (L1,L4)
        dN.row(7) <<   4*L3,  4*L2,  0;        // dN8 (L2,L3)
        dN.row(8) <<   4*L4,   0,   4*L2;      // dN9 (L2,L4)
        dN.row(9) <<   0,    4*L4,  4*L3;      // dN10(L3,L4)
        return dN;
    }
    throw std::invalid_argument("Tetrahedron shape function derivative order > 2 not yet implemented.");
}


} // namespace Utils