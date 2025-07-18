#include "utils/Quadrature.hpp"
#include <stdexcept>
#include <vector> // Ensure vector is included

namespace Utils {

// Note: Points are provided in natural coordinates.
// For a line: [-1, 1]
// For a triangle: (0,0), (1,0), (0,1)
// For a tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)

std::vector<QuadraturePoint> Quadrature::getLineQuadrature(int order) {
    std::vector<QuadraturePoint> points;
    switch (order) {
        case 1: // 1-point, exact for degree 1
            points.push_back({(Eigen::Vector<double, 1>() << 0.0).finished(), 2.0});
            return points;
        case 2: // 2-point, exact for degree 3
            points.push_back({(Eigen::Vector<double, 1>() << -1.0/sqrt(3.0)).finished(), 1.0});
            points.push_back({(Eigen::Vector<double, 1>() <<  1.0/sqrt(3.0)).finished(), 1.0});
            return points;
        case 3: // 3-point, exact for degree 5
            points.push_back({(Eigen::Vector<double, 1>() << -sqrt(3.0/5.0)).finished(), 5.0/9.0});
            points.push_back({(Eigen::Vector<double, 1>() << 0.0).finished(),             8.0/9.0});
            points.push_back({(Eigen::Vector<double, 1>() <<  sqrt(3.0/5.0)).finished(), 5.0/9.0});
            return points;
        case 4: // 4-point, exact for degree 7
            points.push_back({(Eigen::Vector<double, 1>() << -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0))).finished(), (18.0 - sqrt(30.0))/36.0});
            points.push_back({(Eigen::Vector<double, 1>() <<  sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0))).finished(), (18.0 - sqrt(30.0))/36.0});
            points.push_back({(Eigen::Vector<double, 1>() << -sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0))).finished(), (18.0 + sqrt(30.0))/36.0});
            points.push_back({(Eigen::Vector<double, 1>() <<  sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0))).finished(), (18.0 + sqrt(30.0))/36.0});
            return points;
        case 5: // 5-point, exact for degree 9
            points.push_back({(Eigen::Vector<double, 1>() << 0.0).finished(), 128.0/225.0});
            points.push_back({(Eigen::Vector<double, 1>() << -1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0))).finished(), (322.0 + 13.0 * sqrt(70.0))/900.0});
            points.push_back({(Eigen::Vector<double, 1>() <<  1.0/3.0 * sqrt(5.0 - 2.0 * sqrt(10.0/7.0))).finished(), (322.0 + 13.0 * sqrt(70.0))/900.0});
            points.push_back({(Eigen::Vector<double, 1>() << -1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0))).finished(), (322.0 - 13.0 * sqrt(70.0))/900.0});
            points.push_back({(Eigen::Vector<double, 1>() <<  1.0/3.0 * sqrt(5.0 + 2.0 * sqrt(10.0/7.0))).finished(), (322.0 - 13.0 * sqrt(70.0))/900.0});
            return points;
        default:
            throw std::runtime_error("Quadrature order not implemented for Line.");
    }
}

std::vector<QuadraturePoint> Quadrature::getTriangleQuadrature(int order) {
    std::vector<QuadraturePoint> points;
    switch (order) {
        case 1:
            points.push_back({(Eigen::Vector2d() << 1.0/3.0, 1.0/3.0).finished(), 0.5});
            return points;
        case 2:
            points.push_back({(Eigen::Vector2d() << 0.5, 0.0).finished(), 1.0/6.0});
            points.push_back({(Eigen::Vector2d() << 0.5, 0.5).finished(), 1.0/6.0});
            points.push_back({(Eigen::Vector2d() << 0.0, 0.5).finished(), 1.0/6.0});
            return points;
        case 3:
            points.push_back({(Eigen::Vector2d() << 1.0/3.0, 1.0/3.0).finished(), -27.0/96.0});
            points.push_back({(Eigen::Vector2d() << 0.6, 0.2).finished(),      25.0/96.0});
            points.push_back({(Eigen::Vector2d() << 0.2, 0.6).finished(),      25.0/96.0});
            points.push_back({(Eigen::Vector2d() << 0.2, 0.2).finished(),      25.0/96.0});
            return points;
        case 4:
        case 5:
            points.push_back({(Eigen::Vector2d() << 1.0/3.0, 1.0/3.0).finished(), 0.225 * 0.5});
            points.push_back({(Eigen::Vector2d() << 0.797426985353087, 0.101286507323456).finished(), 0.125939180544827 * 0.5});
            points.push_back({(Eigen::Vector2d() << 0.101286507323456, 0.797426985353087).finished(), 0.125939180544827 * 0.5});
            points.push_back({(Eigen::Vector2d() << 0.101286507323456, 0.101286507323456).finished(), 0.125939180544827 * 0.5});
            points.push_back({(Eigen::Vector2d() << 0.059715871789770, 0.470142064105115).finished(), 0.132394152788506 * 0.5});
            points.push_back({(Eigen::Vector2d() << 0.470142064105115, 0.059715871789770).finished(), 0.132394152788506 * 0.5});
            points.push_back({(Eigen::Vector2d() << 0.470142064105115, 0.470142064105115).finished(), 0.132394152788506 * 0.5});
            return points;
        default:
            throw std::runtime_error("Quadrature order not implemented for Triangle.");
    }
}

std::vector<QuadraturePoint> Quadrature::getTetrahedronQuadrature(int order) {
    std::vector<QuadraturePoint> points;
    switch (order) {
        case 1:
            points.push_back({(Eigen::Vector3d() << 0.25, 0.25, 0.25).finished(), 1.0/6.0});
            return points;
        case 2:
            {
                double a = 0.5854101966249685;
                double b = 0.1381966011250105;
                points.push_back({(Eigen::Vector3d() << a, b, b).finished(), 1.0/24.0});
                points.push_back({(Eigen::Vector3d() << b, a, b).finished(), 1.0/24.0});
                points.push_back({(Eigen::Vector3d() << b, b, a).finished(), 1.0/24.0});
                points.push_back({(Eigen::Vector3d() << b, b, b).finished(), 1.0/24.0});
                return points;
            }
        case 3:
        case 4:
        case 5:
            points.push_back({(Eigen::Vector3d() << 0.25, 0.25, 0.25).finished(), -0.8/6.0});
            points.push_back({(Eigen::Vector3d() << 0.5, 1.0/6.0, 1.0/6.0).finished(), 0.45/6.0});
            points.push_back({(Eigen::Vector3d() << 1.0/6.0, 0.5, 1.0/6.0).finished(), 0.45/6.0});
            points.push_back({(Eigen::Vector3d() << 1.0/6.0, 1.0/6.0, 0.5).finished(), 0.45/6.0});
            points.push_back({(Eigen::Vector3d() << 1.0/6.0, 1.0/6.0, 1.0/6.0).finished(), 0.45/6.0});
            return points;
        default:
            throw std::runtime_error("Quadrature order not implemented for Tetrahedron.");
    }
}

} // namespace Utils