#ifndef FEVALUES_HPP
#define FEVALUES_HPP

#include <vector>
#include <Eigen/Dense>
#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"
#include "core/mesh/ElementGeometry.hpp"

namespace Core {

    /**
     * @class FEValues
     * @brief Finite Element Values calculator.
     *
     * For a given element geometry and approximation order, this class pre-calculates
     * and caches all values needed for FEM assembly at all quadrature points. This
     * includes shape function values (N), their gradients in real coordinates (∇N),
     * and the Jacobian determinant multiplied by the quadrature weight (detJ * w_q).
     */
    class FEValues {
    public:
        FEValues(const ElementGeometry& geom, int order, int quad_order);

        // Re-initializes the object to provide values for a specific quadrature point.
        void reinit(int q_point_index);

        // Returns the number of quadrature points.
        size_t num_quadrature_points() const { return quadrature_points_.size(); }

        // --- Accessor Methods for the current quadrature point ---

        // Get shape function values, N
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

        // Pre-calculated values for ALL quadrature points
        const ElementGeometry& geometry_;
        int fe_order_;
        std::vector<Utils::QuadraturePoint> quadrature_points_;
        std::vector<Eigen::VectorXd> all_N_values_;
        std::vector<Eigen::MatrixXd> all_dN_dx_values_;
        std::vector<double> all_detJ_x_w_;
    };

} // namespace Core

#endif // FEVALUES_HPP