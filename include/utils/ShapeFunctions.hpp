#ifndef SHAPEFUNCTIONS_HPP
#define SHAPEFUNCTIONS_HPP
#define MAX_ELEMENT_ORDER_SUPPORTED 2
#include <Eigen/Dense>
#include <vector>

namespace Utils {
    constexpr int MAX_LINE_ORDER_SUPPORTED = 5;
    constexpr int MAX_TRI_ORDER_SUPPORTED = 2;
    constexpr int MAX_TET_ORDER_SUPPORTED = 2;
    class ShapeFunctions {
    public:
        // For 1D Line Elements
        static Eigen::VectorXd getLineShapeFunctions(int order, double xi);
        static Eigen::VectorXd getLineShapeFunctionDerivatives(int order, double xi);

        // For 2D Triangle Elements
        static Eigen::VectorXd getTriShapeFunctions(int order, double xi, double eta);
        static Eigen::MatrixXd getTriShapeFunctionDerivatives(int order, double xi, double eta);

        // For 3D Tetrahedron Elements
        static Eigen::VectorXd getTetShapeFunctions(int order, double xi, double eta, double zeta);
        static Eigen::MatrixXd getTetShapeFunctionDerivatives(int order, double xi, double eta, double zeta);
    };

} // namespace Utils

#endif // SHAPEFUNCTIONS_HPP