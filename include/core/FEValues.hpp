#ifndef FEVALUES_HPP
#define FEVALUES_HPP

#include "mesh/Element.hpp" // 引用 Element 而不是 ElementGeometry
#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"
#include <Eigen/Dense>

namespace Core {

    // Forward declare to break circular dependency
    class Element;

    class FEValues {
    public:
        // FIX: The constructor now takes a const Element&
        FEValues(const Element& elem, int order, int quad_order);

        void reinit(int q_point_index);

        const Eigen::VectorXd& get_shape_values() const;
        const Eigen::MatrixXd& get_shape_gradients() const;
        double get_detJ_times_weight() const;

        size_t num_quadrature_points() const { return quadrature_points_.size(); }

    private:
        // Keep a reference to the element itself
        const Element& element_;
        int element_order_;
        int quadrature_order_;

        std::vector<Utils::QuadraturePoint> quadrature_points_;

        // Values at the current quadrature point
        Eigen::VectorXd N_;
        Eigen::MatrixXd dNdx_;
        double detJ_x_w_;

        // Pre-calculated values for all quadrature points
        std::vector<Eigen::VectorXd> N_values_;
        std::vector<Eigen::MatrixXd> dNdx_values_;
        std::vector<double> detJ_x_w_values_;
    };

} // namespace Core

#endif // FEVALUES_HPP