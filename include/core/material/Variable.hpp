#ifndef VARIABLE_HPP
#define VARIABLE_HPP

#include <string>

namespace Core {

/**
 * @class Variable
 * @brief Represents a field variable in the simulation (e.g., temperature, voltage, pressure)
 *
 * This class encapsulates a single variable that can be used throughout the simulation.
 * Variables are identified by unique names and can have default values.
 */
class Variable {
public:
    /**
     * @brief Constructor
     * @param name Unique name of the variable (e.g., "T" for temperature, "V" for voltage)
     * @param default_value Default value for the variable
     * @param description Human-readable description of the variable
     */
    Variable(const std::string& name, double default_value = 0.0,
             const std::string& description = "");

    /**
     * @brief Get the variable name
     * @return The unique name of the variable
     */
    const std::string& getName() const { return name_; }

    /**
     * @brief Get the default value
     * @return The default value of the variable
     */
    double getDefaultValue() const { return default_value_; }

    /**
     * @brief Get the description
     * @return Human-readable description of the variable
     */
    const std::string& getDescription() const { return description_; }

    /**
     * @brief Set the default value
     * @param value New default value
     */
    void setDefaultValue(double value) { default_value_ = value; }

private:
    std::string name_;
    double default_value_;
    std::string description_;
};

} // namespace Core

#endif // VARIABLE_HPP
