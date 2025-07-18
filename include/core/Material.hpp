#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <string>
#include <map>
#include <any>
#include <functional> // Required for std::function
#include "utils/SimpleLogger.hpp"

namespace Core {

    /**
     * @class Material
     * @brief Represents a material with properties that can be constant or defined by arbitrary user functions.
     */
    class Material {
    public:
        // Type alias for a function that calculates a property based on a set of field values.
        // The map key is the variable name (e.g., "Temperature"), and the value is its magnitude.
        using PropertyFunction = std::function<double(const std::map<std::string, double>&)>;

        explicit Material(const std::string& name);

        // Set a constant scalar property
        void setProperty(const std::string& prop_name, double value);

        // Set a field-dependent property using a function (e.g., a lambda)
        void setProperty(const std::string& prop_name, PropertyFunction func);

        // Get a property's value. This single method handles both constant and function-based properties.
        double getProperty(const std::string& prop_name, const std::map<std::string, double>& field_values = {}) const;

        const std::string& getName() const;

    private:
        std::string name_;
        // std::any allows storing different types: double for constants, PropertyFunction for models.
        std::map<std::string, std::any> properties_;
    };

} // namespace Core

#endif // MATERIAL_HPP