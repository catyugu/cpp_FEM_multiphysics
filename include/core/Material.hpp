#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <string>
#include <map>
#include <any>
#include "utils/SimpleLogger.hpp"

namespace Core {

    /**
     * @class Material
     * @brief Represents a material with properties that can be constant or temperature-dependent.
     */
    class Material {
    public:
        explicit Material(const std::string& name);

        // Set a constant scalar property (e.g., thermal_conductivity)
        void setProperty(const std::string& prop_name, double value);

        // Get a property value. Overloaded for constant and temperature-dependent cases.
        double getProperty(const std::string& prop_name) const;
        double getProperty(const std::string& prop_name, double temperature) const;

        const std::string& getName() const;

    private:
        std::string name_;
        // std::any allows us to store different types of data (double for constant, map for models)
        std::map<std::string, std::any> properties_;
    };

} // namespace Core

#endif // MATERIAL_HPP
