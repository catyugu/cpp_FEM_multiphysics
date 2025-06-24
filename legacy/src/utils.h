// Created by HUAWEI on 2025/6/21.
// Corrected by Gemini
//

#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <stdexcept> // Required for std::logic_error
#include <Eigen/Dense>

// Structure to hold Gaussian Quadrature nodes and weights
// FIX: Corrected the nodes and weights for 20 points and ensured consistency.
// The original implementation had mismatched data which was the primary source of inaccurate results.
// Structure to hold Gaussian Quadrature nodes and weights
// FIX: The nodes and weights have been replaced with a canonical, high-precision
// set of values for 20-point Gauss-Legendre quadrature. The previous set was
// incorrect and caused significant accuracy errors in integration tests.
struct GaussLegendreQuadrature {
    const int num_points;
    const std::vector<double> nodes;
    const std::vector<double> weights;

    GaussLegendreQuadrature() :
        num_points(20),
        nodes({
            -0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188,
            -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154195,
            -0.2277858511416451, -0.0765265211334973,  0.0765265211334973,  0.2277858511416451,
             0.3737060887154195,  0.5108670019508271,  0.6360536807265150,  0.7463319064601508,
             0.8391169718222188,  0.9122344282513259,  0.9639719272779138,  0.9931285991850949
        }),
        weights({
            0.0176140071391521, 0.0406014298003869, 0.0626720483341091, 0.0832767415767047,
            0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183821,
            0.1491729864726037, 0.1527533871307257, 0.1527533871307257, 0.1491729864726037,
            0.1420961093183821, 0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
            0.0832767415767047, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521
        }) {
        if (nodes.size() != num_points || weights.size() != num_points) {
            throw std::logic_error("GaussLegendreQuadrature data mismatch.");
        }
    }
};

const GaussLegendreQuadrature quadrature;

template<typename Func>
double numerical_integral_1D(Func func, double x_min, double x_max) {
    double integral = 0.0;
    const double x_half_diff = (x_max - x_min) / 2.0;
    const double x_half_sum = (x_max + x_min) / 2.0;

    for (int i = 0; i < quadrature.num_points; ++i) {
        integral += quadrature.weights[i] * func(x_half_diff * quadrature.nodes[i] + x_half_sum);
    }

    return integral * x_half_diff;
}

template<typename Func>
double numerical_integral_2D(Func func, double x_min, double x_max, double y_min, double y_max) {
    auto inner_integral = [&](double x) {
        auto func_1d = [&](double y) {
            return func(Eigen::Vector2d(x, y));
        };
        return numerical_integral_1D(func_1d, y_min, y_max);
    };

    return numerical_integral_1D(inner_integral, x_min, x_max);
}

template<typename Func>
double numerical_integral_3D(Func func, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) {
    auto inner_integral = [&](double x, double y) {
        auto func_1d = [&](double z) {
            return func(Eigen::Vector3d(x, y, z));
        };
        return numerical_integral_1D(func_1d, z_min, z_max);
    };

    auto middle_integral = [&](double x) {
        return numerical_integral_1D([&](double y){ return inner_integral(x, y); }, y_min, y_max);
    };
    return numerical_integral_1D(middle_integral, x_min, x_max);
}

/**
 * @brief Computes the line integral of a vector field F along a parameterized path r(t).
 * The integral is evaluated as: integral from t_min to t_max of F(r(t)) . r'(t) dt.
 *
 * @tparam Func1 A callable that takes an Eigen::Vector3d (position) and returns an Eigen::Vector3d (field vector).
 * @tparam Func2 A callable that takes a double (t) and returns an Eigen::Vector3d (position vector).
 * @param func The vector field F.
 * @param path The parameterized path r(t).
 * @param t_min The starting parameter for the path.
 * @param t_max The ending parameter for the path.
 * @return The approximate value of the line integral.
 */
template<typename Func1, typename Func2>
double numerical_integral_path(Func1 func, Func2 path, double t_min, double t_max) {
    // A small step for numerical differentiation.
    // NOTE: This value may need tuning for functions with very high or low frequency changes.
    const double h = 1e-8;

    auto integrand = [&](double t) {
        // Numerical derivative of the path using a central difference scheme for better accuracy.
        Eigen::Vector3d path_derivative = (path(t + h) - path(t - h)) / (2.0 * h);

        return func(path(t)).dot(path_derivative);
    };

    return numerical_integral_1D(integrand, t_min, t_max);
}

#endif //UTILS_H