#ifndef INTERPOLATION_UTILITIES_HPP
#define INTERPOLATION_UTILITIES_HPP

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>
#include "physics/PhysicsField.hpp"

namespace Core {
    class Element;

    class Material;
}

namespace Utils {

/**
 * @class InterpolationUtilities
 * @brief Utilities for interpolating field values at quadrature points
 *
 * This class provides methods to interpolate field variables at specific
 * quadrature points using shape functions, enabling node-wise material
 * property calculations as used in commercial software.
 */
class InterpolationUtilities {
public:

    /**
     * @brief Simple interpolation of variables at a quadrature point using shape functions
     * @param element Element to interpolate within
     * @param shape_values Shape function values at the quadrature point
     * @param variable_names List of variables to interpolate
     * @param physics_fields Map of physics fields to get solution values from
     * @return Map of variable names to interpolated values at the quadrature point
     */
    static std::map<std::string, double> interpolateAtQuadraturePoint(
        const Core::Element* element,
        const Eigen::VectorXd& shape_values,
        const std::vector<std::string>& variable_names,
        const std::map<std::string, const Physics::PhysicsField*>& physics_fields
    );

    /**
     * @brief Evaluate material property at a quadrature point using interpolated variables
     * @param material Material to evaluate
     * @param property_name Name of the property to evaluate
     * @param interpolated_variables Map of interpolated variable values
     * @return Evaluated material property value
     */
    static double evaluateMaterialPropertyAtQuadraturePoint(
        const Core::Material& material,
        const std::string& property_name,
        const std::map<std::string, double>& interpolated_variables
    );

private:
    InterpolationUtilities() = default;
};

} // namespace Utils

#endif // INTERPOLATION_UTILITIES_HPP
