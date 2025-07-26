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

        if (num_vertices == 2) { // P2 Line: v0, v1, midpoint
            all_coords.col(0) = vertex_coords.col(0);
            all_coords.col(1) = vertex_coords.col(1);
            all_coords.col(2) = (vertex_coords.col(0) + vertex_coords.col(1)) * 0.5;
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


// 新的、轻量化的构造函数
FEValues::FEValues(const ElementGeometry& geom, int order, const ReferenceElementData& ref_data)
    : geometry_(geom), fe_order_(order), ref_data_(ref_data) {

    Eigen::MatrixXd all_node_coords = calculate_all_node_coords(geometry_, fe_order_);

    // 直接使用缓存的参考数据进行几何映射计算
    for (size_t qp_idx = 0; qp_idx < ref_data_.quadrature_points.size(); ++qp_idx) {
        const auto& dN_d_natural = ref_data_.dN_d_natural_at_qps[qp_idx];

        Eigen::MatrixXd J = all_node_coords * dN_d_natural;
        double detJ = J.determinant();

        if (detJ <= 1e-20) {
            throw Exception::SolverException("Zero or negative Jacobian determinant encountered.");
        }
        Eigen::MatrixXd J_inv = J.inverse();
        Eigen::MatrixXd dN_dx = dN_d_natural * J_inv;

        all_dN_dx_values_.push_back(dN_dx.transpose());
        all_detJ_x_w_.push_back(detJ * ref_data_.quadrature_points[qp_idx].weight);
    }

    if (!ref_data_.quadrature_points.empty()) {
        reinit(0);
    }
}

void FEValues::reinit(int q_point_index) {
    if (q_point_index >= num_quadrature_points()) {
        throw std::out_of_range("Quadrature point index is out of range.");
    }
    // N_values_ 直接从缓存的参考数据中获取
    N_values_ = ref_data_.N_values_at_qps[q_point_index];
    dN_dx_values_ = all_dN_dx_values_[q_point_index];
    detJ_x_w_ = all_detJ_x_w_[q_point_index];
}

} // namespace Core