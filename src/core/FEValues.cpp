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

        // --- THIS IS THE DEFINITIVE FIX ---
        // Assemble all coordinate vectors according to the strict canonical ordering.
        if (num_vertices == 2) { // P2 Line: v0, midpoint, v1
            all_coords.col(0) = vertex_coords.col(0);
            all_coords.col(1) = (vertex_coords.col(0) + vertex_coords.col(1)) * 0.5;
            all_coords.col(2) = vertex_coords.col(1);
        } else {
            // P2 Tri/Tet: Vertices first, then edge midpoints in canonical order.
            all_coords.leftCols(num_vertices) = vertex_coords;
            int edge_node_idx = num_vertices;
            std::vector<std::pair<int, int>> edges;
            if (num_vertices == 3)      edges = {{0,1}, {1,2}, {2,0}};                      // Tri6
            else if (num_vertices == 4) edges = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}}; // Tet10

            for (const auto& edge : edges) {
                all_coords.col(edge_node_idx++) = (vertex_coords.col(edge.first) + vertex_coords.col(edge.second)) * 0.5;
            }
        }
        return all_coords;
    }

    throw std::runtime_error("FEValues currently only supports order 1 and 2.");
}


FEValues::FEValues(const ElementGeometry& geom, int order, int quad_order)
    : geometry_(geom), fe_order_(order) {

    Eigen::MatrixXd all_node_coords = calculate_all_node_coords(geometry_, fe_order_);

    switch(geom.get_num_vertices()) {
        case 2: quadrature_points_ = Utils::Quadrature::getLineQuadrature(quad_order); break;
        case 3: quadrature_points_ = Utils::Quadrature::getTriangleQuadrature(quad_order); break;
        case 4: quadrature_points_ = Utils::Quadrature::getTetrahedronQuadrature(quad_order); break;
        default: throw Exception::ConfigurationException("Unsupported element type for FEValues.");
    }

    for (const auto& qp : quadrature_points_) {
        Eigen::VectorXd N_natural;
        Eigen::MatrixXd dN_d_natural;

        switch(geom.get_num_vertices()) {
            case 2:
                N_natural = Utils::ShapeFunctions::getLineShapeFunctions(fe_order_, qp.point(0));
                dN_d_natural = Utils::ShapeFunctions::getLineShapeFunctionDerivatives(fe_order_, qp.point(0));
                break;
            case 3:
                N_natural = Utils::ShapeFunctions::getTriShapeFunctions(fe_order_, qp.point(0), qp.point(1));
                dN_d_natural = Utils::ShapeFunctions::getTriShapeFunctionDerivatives(fe_order_, qp.point(0), qp.point(1));
                break;
            case 4:
                N_natural = Utils::ShapeFunctions::getTetShapeFunctions(fe_order_, qp.point(0), qp.point(1), qp.point(2));
                dN_d_natural = Utils::ShapeFunctions::getTetShapeFunctionDerivatives(fe_order_, qp.point(0), qp.point(1), qp.point(2));
                break;
        }

        Eigen::MatrixXd J = all_node_coords * dN_d_natural;
        double detJ = J.determinant();

        if (detJ <= 1e-20) {
            throw Exception::SolverException("Zero or negative Jacobian determinant encountered.");
        }
        Eigen::MatrixXd J_inv = J.inverse();
        Eigen::MatrixXd dN_dx = dN_d_natural * J_inv;

        all_N_values_.push_back(N_natural);
        all_dN_dx_values_.push_back(dN_dx.transpose());
        all_detJ_x_w_.push_back(detJ * qp.weight);
    }

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