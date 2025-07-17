#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"
#include <cmath>

namespace Core {

// Helper function to calculate the coordinates of all nodes (vertices + higher-order)
// This is crucial for correct isoparametric mapping in higher-order elements.
static Eigen::MatrixXd calculate_all_node_coords(const Core::ElementGeometry& geom, int order) {
    const auto& vertex_coords = geom.get_vertex_coords();
    const size_t num_vertices = geom.get_num_vertices();
    const int dim = geom.get_dimension();

    if (order == 1) {
        return vertex_coords;
    }

    if (order == 2) {
        size_t num_elem_nodes;
        if (num_vertices == 2) num_elem_nodes = 3;      // Line3
        else if (num_vertices == 3) num_elem_nodes = 6; // Tri6
        else if (num_vertices == 4) num_elem_nodes = 10;// Tet10
        else throw std::runtime_error("Unsupported element type for order 2 in FEValues coord calculation.");

        Eigen::MatrixXd all_coords(dim, num_elem_nodes);
        // First, copy the vertex coordinates
        all_coords.leftCols(num_vertices) = vertex_coords;

        // Next, calculate the midpoint coordinates for the edges
        int edge_node_idx = num_vertices;
        std::vector<std::pair<int, int>> edges;
        if (num_vertices == 2) edges = {{0,1}}; // Line
        if (num_vertices == 3) edges = {{0,1}, {1,2}, {2,0}}; // Triangle
        if (num_vertices == 4) edges = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}}; // Tetrahedron

        for (const auto& edge : edges) {
            all_coords.col(edge_node_idx++) = (vertex_coords.col(edge.first) + vertex_coords.col(edge.second)) * 0.5;
        }
        return all_coords;
    }

    throw std::runtime_error("FEValues currently only supports order 1 and 2.");
}


FEValues::FEValues(const ElementGeometry& geom, int order, int quad_order)
    : geometry_(geom), fe_order_(order) {

    // 1. Calculate coordinates of ALL nodes required for the given order.
    Eigen::MatrixXd all_node_coords = calculate_all_node_coords(geometry_, fe_order_);

    // 2. Get the correct quadrature rule based on the element type and order.
    switch(geom.get_num_vertices()) {
        case 2: // Line
            quadrature_points_ = Utils::Quadrature::getLineQuadrature(quad_order);
            break;
        case 3: // Triangle
            quadrature_points_ = Utils::Quadrature::getTriangleQuadrature(quad_order);
            break;
        case 4: // Tetrahedron
            quadrature_points_ = Utils::Quadrature::getTetrahedronQuadrature(quad_order);
            break;
        default:
            throw Exception::ConfigurationException("Unsupported element type for FEValues.");
    }

    // 3. Pre-calculate and cache all values for all quadrature points.
    for (const auto& qp : quadrature_points_) {
        // Get shape functions and their derivatives in NATURAL coordinates
        Eigen::VectorXd N_natural;
        Eigen::MatrixXd dN_d_natural;

        switch(geom.get_num_vertices()) {
            case 2: // Line
                N_natural = Utils::ShapeFunctions::getLineShapeFunctions(fe_order_, qp.point(0));
                dN_d_natural = Utils::ShapeFunctions::getLineShapeFunctionDerivatives(fe_order_, qp.point(0));
                break;
            case 3: // Triangle
                N_natural = Utils::ShapeFunctions::getTriShapeFunctions(fe_order_, qp.point(0), qp.point(1));
                dN_d_natural = Utils::ShapeFunctions::getTriShapeFunctionDerivatives(fe_order_, qp.point(0), qp.point(1));
                break;
            case 4: // Tetrahedron
                N_natural = Utils::ShapeFunctions::getTetShapeFunctions(fe_order_, qp.point(0), qp.point(1), qp.point(2));
                dN_d_natural = Utils::ShapeFunctions::getTetShapeFunctionDerivatives(fe_order_, qp.point(0), qp.point(1), qp.point(2));
                break;
        }

        // Calculate Jacobian using ALL node coordinates (vertices + higher-order)
        // J_ij = sum_k (X_ik * dNk/d_xi_j) => J = all_node_coords * dN_d_natural
        Eigen::MatrixXd J = all_node_coords * dN_d_natural;
        double detJ = J.determinant(); // This will now operate on a square matrix.

        if (detJ <= 1e-12) {
            throw Exception::SolverException("Zero or negative Jacobian determinant encountered.");
        }
        Eigen::MatrixXd J_inv = J.inverse();

        // Calculate shape gradients in REAL coordinates (∇N = J^-T * ∇_nat(N))
        // or B = (J^-1 * (dN/d_xi)^T)^T = (dN/d_xi) * J_inv
        Eigen::MatrixXd dN_dx = dN_d_natural * J_inv;

        // Cache the results for this quadrature point
        all_N_values_.push_back(N_natural);
        all_dN_dx_values_.push_back(dN_dx.transpose()); // Store as B-matrix format
        all_detJ_x_w_.push_back(detJ * qp.weight);
    }

    // Initialize the current values to the first quadrature point
    if (!quadrature_points_.empty()) {
        reinit(0);
    }
}

void FEValues::reinit(int q_point_index) {
    if (q_point_index >= num_quadrature_points()) {
        throw std::out_of_range("Quadrature point index is out of range.");
    }
    N_values_ = all_N_values_[q_point_index];
    dN_dx_values_ = all_dN_dx_values_[q_point_index];
    detJ_x_w_ = all_detJ_x_w_[q_point_index];
}

} // namespace Core