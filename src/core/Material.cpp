#include "core/Material.hpp"
#include <stdexcept>

namespace Core {

Material::Material(const std::string& name) : name_(name) {}

void Material::setProperty(const std::string& prop_name, double value) {
    properties_[prop_name] = value;
}

const std::string& Material::getName() const {
    return name_;
}

// Overload for constant properties
double Material::getProperty(const std::string& prop_name) const {
    auto it = properties_.find(prop_name);
    if (it == properties_.end()) {
        throw std::runtime_error("Material property '" + prop_name + "' not found.");
    }
    try {
        return std::any_cast<double>(it->second);
    } catch (const std::bad_any_cast& e) {
        throw std::runtime_error("Property '" + prop_name + "' is not a constant value.");
    }
}

// Overload for temperature-dependent properties
double Material::getProperty(const std::string& prop_name, double temperature) const {
    auto it = properties_.find(prop_name);
    if (it == properties_.end()) {
        throw std::runtime_error("Material property '" + prop_name + "' not found.");
    }

    // Try to cast to a map for model parameters
    try {
        const auto& params = std::any_cast<const std::map<std::string, double>&>(it->second);

        // --- Linear Temperature Model for Electrical Conductivity ---
        // This is where you would implement different material models.
        if (prop_name == "electrical_conductivity") {
            double sigma_ref = params.at("sigma_ref");
            double alpha = params.at("alpha");
            double T_ref = params.at("T_ref");
            return sigma_ref * (1 + alpha * (temperature - T_ref));
        }

        // Add other models here...

        throw std::runtime_error("Unknown temperature-dependent model for property: " + prop_name);

    } catch (const std::bad_any_cast&) {
        // If it's not a map, it must be a constant property, so we return that.
        // This allows thermal conductivity (k) to be defined as a constant
        // while electrical conductivity (sigma) is temperature-dependent.
        return std::any_cast<double>(it->second);
    } catch (const std::out_of_range& e) {
        throw std::runtime_error("Missing parameter in model for property '" + prop_name + "'.");
    }
}

} // namespace Core
