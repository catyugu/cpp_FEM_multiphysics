#ifndef SHAPEFUNCTIONS_HPP
#define SHAPEFUNCTIONS_HPP

#include <Eigen/Dense>
#include <vector>

namespace Utils {

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