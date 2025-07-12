#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>
#include <Eigen/Dense>

namespace Utils {

    struct QuadraturePoint {
        Eigen::VectorXd point;
        double weight;
    };

    class Quadrature {
    public:
        static std::vector<QuadraturePoint> getLineQuadrature(int order);
        static std::vector<QuadraturePoint> getTriangleQuadrature(int order);
        static std::vector<QuadraturePoint> getTetrahedronQuadrature(int order);
    };

} // namespace Utils

#endif // QUADRATURE_HPP