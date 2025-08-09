#ifndef VARIABLE_MANAGER_HPP
#define VARIABLE_MANAGER_HPP

#include "Variable.hpp"
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <stdexcept>

namespace Core {

/**
 * @class VariableManager
 * @brief Global manager for all simulation variables
 *
 * This singleton class manages all variables in the simulation, providing
 * registration, retrieval, and value management functionality. It serves as
 * the central registry for all field variables used in multiphysics simulations.
 */
class VariableManager {
public:
    /**
     * @brief Get the singleton instance
     * @return Reference to the global VariableManager instance
     */
    static VariableManager& getInstance();

    /**
     * @brief Register a new variable
     * @param name Unique name for the variable
     * @param default_value Default value for the variable
     * @param description Human-readable description
     * @return Reference to the created Variable
     * @throw std::invalid_argument if variable name already exists
     */
    const Variable& registerVariable(const std::string& name,
                                   double default_value = 0.0,
                                   const std::string& description = "");

    /**
     * @brief Get a variable by name
     * @param name Name of the variable to retrieve
     * @return Reference to the Variable
     * @throw std::out_of_range if variable not found
     */
    const Variable& getVariable(const std::string& name) const;

    /**
     * @brief Check if a variable exists
     * @param name Name of the variable to check
     * @return True if variable exists, false otherwise
     */
    bool hasVariable(const std::string& name) const;

    /**
     * @brief Get all registered variables
     * @return Vector of all variable names
     */
    std::vector<std::string> getAllVariableNames() const;

    /**
     * @brief Clear all registered variables (mainly for testing)
     */
    void clear();

    /**
     * @brief Get the number of registered variables
     * @return Number of variables
     */
    size_t getVariableCount() const { return variables_.size(); }

private:
    // Private constructor for singleton
    VariableManager() = default;

    // Delete copy constructor and assignment operator
    VariableManager(const VariableManager&) = delete;
    VariableManager& operator=(const VariableManager&) = delete;

    std::map<std::string, std::unique_ptr<Variable>> variables_;
};

} // namespace Core

#endif // VARIABLE_MANAGER_HPP
