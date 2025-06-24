#include "core/Material.hpp"
#include <stdexcept>

namespace Core {

    Material::Material(const std::string& name) : name_(name) {}

    void Material::setProperty(const std::string& prop_name, double value) {
        properties_[prop_name] = value;
    }

    double Material::getProperty(const std::string& prop_name) const {
        auto it = properties_.find(prop_name);
        if (it == properties_.end()) {
            auto msg = "Material '" + name_ + "' does not have property '" + prop_name + "'.";
            SimpleLogger::Logger::instance().error(msg);
            throw std::runtime_error(msg);
        }
        return it->second;
    }

    bool Material::hasProperty(const std::string& prop_name) const {
        return properties_.find(prop_name) != properties_.end();
    }

    const std::string& Material::getName() const {
        return name_;
    }

} // namespace Core
