#include "utils/Quadrature.hpp"
#include <stdexcept>

namespace Utils {

    std::vector<QuadraturePoint> Quadrature::getLineQuadrature(int order) {
        if (order == 1) {
            return {{Eigen::Vector<double, 1>(0.0), 2.0}};
        }
        throw std::runtime_error("Quadrature order not implemented for Line.");
    }

    std::vector<QuadraturePoint> Quadrature::getTriangleQuadrature(int order) {
        if (order == 1) {
            return {{Eigen::Vector2d(1.0/3.0, 1.0/3.0), 0.5}};
        }
        throw std::runtime_error("Quadrature order not implemented for Triangle.");
    }

    std::vector<QuadraturePoint> Quadrature::getTetrahedronQuadrature(int order) {
        if (order == 1) {
            return {{Eigen::Vector3d(0.25, 0.25, 0.25), 1.0/6.0}};
        }
        throw std::runtime_error("Quadrature order not implemented for Tetrahedron.");
    }

} // namespace Utils