#ifndef FEVALUES_HPP
#define FEVALUES_HPP

#include <vector>
#include <Eigen/Dense>
#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"
#include "core/mesh/ElementGeometry.hpp"
#include "core/ReferenceElement.hpp" // 引入新的头文件

namespace Core {

    /**
     * @class FEValues
     * @brief 有限元值计算器 (已重构为轻量级单元视图).
     *
     * 该类接收预计算的参考单元数据和具体单元的几何信息，
     * 实时计算与该单元相关的几何映射值，如真实坐标下的梯度 (∇N)
     * 和雅可比行列式乘积分权重 (detJ * w_q)。
     */
    class FEValues {
    public:
        // 新的构造函数，接收几何信息和缓存的参考单元数据
        FEValues(const ElementGeometry& geom, int order, const ReferenceElementData& ref_data);

        // Re-initializes the object to provide values for a specific quadrature point.
        void reinit(int q_point_index);

        // Returns the number of quadrature points.
        size_t num_quadrature_points() const { return ref_data_.quadrature_points.size(); }

        // --- Accessor Methods for the current quadrature point ---

        // Get shape function values, N (直接从缓存获取)
        const Eigen::VectorXd& get_shape_values() const { return N_values_; }

        // Get shape function gradients in the REAL coordinate system, ∇N
        const Eigen::MatrixXd& get_shape_gradients() const { return dN_dx_values_; }

        // Get the combined Jacobian determinant and quadrature weight, det(J) * w_q
        double get_detJ_times_weight() const { return detJ_x_w_; }

    private:
        // Cached values for the current quadrature point
        Eigen::VectorXd N_values_;
        Eigen::MatrixXd dN_dx_values_;
        double detJ_x_w_;

        // Pre-calculated values for ALL quadrature points of THIS specific element
        const ElementGeometry& geometry_;
        int fe_order_;
        const ReferenceElementData& ref_data_; // 指向缓存数据的引用
        std::vector<Eigen::MatrixXd> all_dN_dx_values_; // 只存储与几何相关的真实坐标导数
        std::vector<double> all_detJ_x_w_;
    };

} // namespace Core

#endif // FEVALUES_HPP