#include "core/Material.hpp"
#include <stdexcept>

namespace Core {

    Material::Material(const std::string& name) : name_(name) {}

    void Material::setProperty(const std::string& prop_name, double value) {
        properties_[prop_name] = value;
    }

    void Material::setProperty(const std::string& prop_name, PropertyFunction func) {
        properties_[prop_name] = func;
    }

    const std::string& Material::getName() const {
        return name_;
    }

    double Material::getProperty(const std::string& prop_name, const std::map<std::string, double>& field_values) const {
        auto it = properties_.find(prop_name);
        if (it == properties_.end()) {
            throw std::runtime_error("Material property '" + prop_name + "' not found.");
        }

        const std::any& prop = it->second;

        // Check if the stored property is a constant double
        if (prop.type() == typeid(double)) {
            return std::any_cast<double>(prop);
        }
        // Check if it is a user-defined function
        else if (prop.type() == typeid(PropertyFunction)) {
            const auto& func = std::any_cast<const PropertyFunction&>(prop);
            return func(field_values); // Execute the function
        }

        throw std::runtime_error("Unknown type for property '" + prop_name + "'.");
    }

} // namespace Core