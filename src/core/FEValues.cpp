#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"

namespace Core {
    static Eigen::MatrixXd calculate_all_node_coords(const Core::ElementGeometry& geom, int order) {
        const auto& vertex_coords = geom.get_vertex_coords();
        const size_t num_vertices = geom.get_num_vertices();
        const int dim = geom.get_dimension();

        if (order == 1) {
            return vertex_coords;
        }

        if (order == 2) {
            size_t num_elem_nodes;
            if (num_vertices == 2) num_elem_nodes = 3;
            else if (num_vertices == 3) num_elem_nodes = 6;
            else if (num_vertices == 4) num_elem_nodes = 10;
            else throw std::runtime_error("Unsupported element type for order 2 in FEValues coord calculation.");

            Eigen::MatrixXd all_coords(dim, num_elem_nodes);

            // Vertices first for all element types
            all_coords.leftCols(static_cast<Eigen::Index>(num_vertices)) = vertex_coords;
            int edge_node_idx = static_cast<int>(num_vertices);

            if (num_vertices == 2) { // P2 Line: v0, v1, midpoint
                all_coords.col(edge_node_idx++) = (vertex_coords.col(0) + vertex_coords.col(1)) * 0.5;
            } else {
                // P2 Tri/Tet: edge midpoints in canonical order.
                std::vector<std::pair<int, int>> edges;
                if (num_vertices == 3) {
                    edges = {{0,1}, {1,2}, {2,0}};
                } else { // num_vertices == 4
                    edges = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
                }

                for (const auto& edge : edges) {
                    all_coords.col(edge_node_idx++) = (vertex_coords.col(edge.first) + vertex_coords.col(edge.second)) * 0.5;
                }
            }
            return all_coords;
        }

        throw std::runtime_error("FEValues currently only supports order 1 and 2.");
    }


    FEValues::FEValues(const ElementGeometry& geom, int order, const ReferenceElementData& ref_data)
        : detJ_x_w_(0.0), geometry_(geom), fe_order_(order), ref_data_(ref_data),
          analysis_type_(AnalysisType::CUSTOM) { // 初始化分析类型为CUSTOM

        Eigen::MatrixXd all_node_coords = calculate_all_node_coords(geometry_, fe_order_);

        for (size_t qp_idx = 0; qp_idx < ref_data_.quadrature_points.size(); ++qp_idx) {
            const auto& dN_d_natural = ref_data_.dN_d_natural_at_qps[qp_idx];

            Eigen::MatrixXd J = all_node_coords * dN_d_natural;
            double detJ = J.determinant();

            if (detJ <= 1e-20) {
                throw Exception::SolverException("Zero or negative Jacobian determinant encountered.");
            }
            Eigen::MatrixXd J_inv = J.inverse();
            Eigen::MatrixXd dN_dx = dN_d_natural * J_inv;

            all_dN_dx_values_.emplace_back(dN_dx.transpose());
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
        N_values_ = ref_data_.N_values_at_qps[q_point_index];
        dN_dx_values_ = all_dN_dx_values_[q_point_index];
        detJ_x_w_ = all_detJ_x_w_[q_point_index];

        // 如果已设置分析类型，更新当前积分点的B矩阵
        if (analysis_type_ != AnalysisType::CUSTOM && q_point_index < all_B_matrices_.size()) {
            B_matrix_ = all_B_matrices_[q_point_index];
        }
    }

    void FEValues::setAnalysisType(AnalysisType type) {
        analysis_type_ = type;

        // 为所有积分点预计算B矩阵
        all_B_matrices_.clear();
        all_B_matrices_.reserve(num_quadrature_points());

        for (size_t qp_idx = 0; qp_idx < num_quadrature_points(); ++qp_idx) {
            buildBMatrix(static_cast<int>(qp_idx));
        }

        // 更新当前积分点的B矩阵
        if (!all_B_matrices_.empty()) {
            B_matrix_ = all_B_matrices_[0];
        }
    }

    void FEValues::buildBMatrix(int q_point_index) {
        const auto& dN_dx = all_dN_dx_values_[q_point_index];

        Eigen::MatrixXd B;
        switch (analysis_type_) {
            case AnalysisType::SCALAR_DIFFUSION:
                B = buildScalarDiffusionBMatrix(dN_dx);
            break;
            case AnalysisType::VECTOR_CURL:
                B = buildVectorCurlBMatrix(dN_dx);
            break;
            case AnalysisType::VECTOR_GRADIENT:
                // TODO: 实现固体力学的应变位移矩阵
                    throw std::runtime_error("VECTOR_GRADIENT B matrix not implemented yet");
            break;
            default:
                throw std::runtime_error("Unknown analysis type for B matrix construction");
        }

        all_B_matrices_.push_back(B);
    }

    Eigen::MatrixXd FEValues::buildScalarDiffusionBMatrix(const Eigen::MatrixXd& dN_dx) {
        // 对于标量扩散问题，B矩阵就是形函数的梯度 ∇N
        // dN_dx 的形状是 (dim, num_nodes)，我们需要转置为 (dim, num_nodes)
        return dN_dx;
    }

    Eigen::MatrixXd FEValues::buildVectorCurlBMatrix(const Eigen::MatrixXd& dN_dx) {
        // 对于3D磁场问题，构建旋度矩阵 B_curl
        // curl(N) = [∂N/∂y*ez - ∂N/∂z*ey, ∂N/∂z*ex - ∂N/∂x*ez, ∂N/∂x*ey - ∂N/∂y*ex]

        const auto num_nodes = dN_dx.cols();
        const auto dim = dN_dx.rows();

        if (dim != 3) {
            throw std::runtime_error("Vector curl B matrix only supported for 3D problems");
        }

        // B_curl 矩阵的大小：(3, num_nodes * 3)
        Eigen::MatrixXd B_curl(3, num_nodes * 3);
        B_curl.setZero();

        for (Eigen::Index i = 0; i < num_nodes; ++i) {
            double dN_dx_val = dN_dx(0, i);  // ∂N_i/∂x
            double dN_dy_val = dN_dx(1, i);  // ∂N_i/∂y
            double dN_dz_val = dN_dx(2, i);  // ∂N_i/∂z

            // curl(N_i * ex) = [0, ∂N_i/∂z, -∂N_i/∂y]
            B_curl(1, i * 3 + 0) = dN_dz_val;
            B_curl(2, i * 3 + 0) = -dN_dy_val;

            // curl(N_i * ey) = [-∂N_i/∂z, 0, ∂N_i/∂x]
            B_curl(0, i * 3 + 1) = -dN_dz_val;
            B_curl(2, i * 3 + 1) = dN_dx_val;

            // curl(N_i * ez) = [∂N_i/∂y, -∂N_i/∂x, 0]
            B_curl(0, i * 3 + 2) = dN_dy_val;
            B_curl(1, i * 3 + 2) = -dN_dx_val;
        }

        return B_curl;
    }
}

