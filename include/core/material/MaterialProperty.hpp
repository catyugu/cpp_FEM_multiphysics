#ifndef MATERIAL_PROPERTY_HPP
#define MATERIAL_PROPERTY_HPP

#include <functional>
#include <string>
#include <map>
#include <vector>
#include <memory>

namespace Core {
    class Element; // Forward declaration

/**
 * @class MaterialProperty
 * @brief Represents a material property that can be constant or depend on field variables
 *
 * This class encapsulates a material property (e.g., thermal conductivity, electrical conductivity)
 * that can be either a constant value or a function of one or more field variables.
 * The property can evaluate itself based on current variable values in an element.
 */
class MaterialProperty {
public:
    /**
     * @brief Function type for variable-dependent properties
     * Takes a map of variable names to values and returns the property value
     */
    using PropertyFunction = std::function<double(const std::map<std::string, double>&)>;

    /**
     * @brief Constructor for constant property
     * @param name Name of the property (e.g., "thermal_conductivity", "electrical_conductivity")
     * @param value Constant value of the property
     */
    MaterialProperty(const std::string& name, double value);

    /**
     * @brief Constructor for variable-dependent property
     * @param name Name of the property
     * @param function Function that computes the property based on variable values
     * @param dependent_variables List of variable names this property depends on
     */
    MaterialProperty(const std::string& name,
                    PropertyFunction function,
                    const std::vector<std::string>& dependent_variables);

    /**
     * @brief Get the property name
     * @return Name of the property
     */
    const std::string& getName() const { return name_; }

    /**
     * @brief Check if this property is constant
     * @return True if constant, false if variable-dependent
     */
    bool isConstant() const { return is_constant_; }

    /**
     * @brief Get the list of variables this property depends on
     * @return Vector of variable names
     */
    const std::vector<std::string>& getDependentVariables() const { return dependent_variables_; }

    /**
     * @brief Evaluate the property for a given element
     * @param element Element to get variable values from
     * @return Computed property value
     */
    double evaluate(const Element* element) const;

    /**
     * @brief Evaluate the property with explicit variable values
     * @param variable_values Map of variable names to values
     * @return Computed property value
     */
    double evaluate(const std::map<std::string, double>& variable_values) const;

    /**
     * @brief Set a constant value (converts to constant property)
     * @param value New constant value
     */
    void setConstantValue(double value);

    /**
     * @brief Set a function (converts to variable-dependent property)
     * @param function New property function
     * @param dependent_variables Variables the function depends on
     */
    void setFunction(PropertyFunction function,
                    const std::vector<std::string>& dependent_variables);

private:
    std::string name_;
    bool is_constant_;
    double constant_value_;
    PropertyFunction function_;
    std::vector<std::string> dependent_variables_;
};

} // namespace Core

#endif // MATERIAL_PROPERTY_HPP
