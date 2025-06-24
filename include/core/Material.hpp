#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <string>
#include <map>
#include "utils/SimpleLogger.hpp"

namespace Core {

    // Represents a material with a set of physical properties.
    class Material {
    public:
        explicit Material(const std::string& name);

        // Set a scalar property (e.g., "thermal_conductivity")
        void setProperty(const std::string& prop_name, double value);

        // Get a property. Throws an exception if the property doesn't exist.
        double getProperty(const std::string& prop_name) const;

        // Check if a property exists
        bool hasProperty(const std::string& prop_name) const;

        const std::string& getName() const;

    private:
        std::string name_;
        std::map<std::string, double> properties_;
    };

} // namespace Core

#endif // MATERIAL_HPP
