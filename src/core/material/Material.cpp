#include "core/material/Material.hpp"
#include "core/material/MaterialProperty.hpp"
#include "core/mesh/Element.hpp"
#include <stdexcept>
#include <utility>

namespace Core {

    Material::Material(int id, std::string  name) : name_(std::move(name)), id_(id) {}

    void Material::setProperty(const std::string& prop_name, double value) {
        properties_[prop_name] = value;
    }

    void Material::setProperty(const std::string& prop_name, PropertyFunction func) {
        properties_[prop_name] = func;
    }

    const std::string& Material::getName() const {
        return name_;
    }



    // ========= 新增：扩展的材料属性管理功能实现 =========

    void Material::setMaterialProperty(const std::string& prop_name,
                                     const MaterialProperty& material_prop) {
        // 将MaterialProperty对象存储到properties_映射中
        properties_[prop_name] = material_prop;
    }

    double Material::getPropertyFromElement(const std::string& prop_name,
                                          const Element* element) const {
        auto it = properties_.find(prop_name);
        if (it == properties_.end()) {
            throw std::runtime_error("Material property '" + prop_name + "' not found.");
        }

        const std::any& prop = it->second;

        // 检查是否为常数double
        if (prop.type() == typeid(double)) {
            return std::any_cast<double>(prop);
        }
        // 检查是否为用户定义的函数
        else if (prop.type() == typeid(PropertyFunction)) {
            if (!element) {
                throw std::invalid_argument("Element cannot be null for function-based property");
            }
            // 从元素获取所有变量值
            const auto& variable_values = element->getAllVariableValues();
            const auto& func = std::any_cast<const PropertyFunction&>(prop);
            return func(variable_values);
        }
        // 检查是否为MaterialProperty对象
        else if (prop.type() == typeid(MaterialProperty)) {
            const auto& material_prop = std::any_cast<const MaterialProperty&>(prop);
            return material_prop.evaluate(element);
        }

        throw std::runtime_error("Unknown type for property '" + prop_name + "'.");
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
        // Check if it is a MaterialProperty object
        else if (prop.type() == typeid(MaterialProperty)) {
            const auto& material_prop = std::any_cast<const MaterialProperty&>(prop);
            return material_prop.evaluate(field_values);
        }

        throw std::runtime_error("Unknown type for property '" + prop_name + "'.");
    }
    double Material::getPropertyAtQuadraturePoint(const std::string& prop_name,
                                                 const std::map<std::string, double>& interpolated_values) const {
        auto it = properties_.find(prop_name);
        if (it == properties_.end()) {
            throw std::runtime_error("Material property '" + prop_name + "' not found.");
        }

        const std::any& prop = it->second;

        // 检查是否为常数double
        if (prop.type() == typeid(double)) {
            return std::any_cast<double>(prop);
        }
        // 检查是否为用户定义的函数
        else if (prop.type() == typeid(PropertyFunction)) {
            const auto& func = std::any_cast<const PropertyFunction&>(prop);
            return func(interpolated_values);
        }
        // 检查是否为MaterialProperty对象
        else if (prop.type() == typeid(MaterialProperty)) {
            const auto& material_prop = std::any_cast<const MaterialProperty&>(prop);
            return material_prop.evaluate(interpolated_values);
        }

        throw std::runtime_error("Unknown type for property '" + prop_name + "'.");
    }

    bool Material::hasProperty(const std::string& prop_name) const {
        return properties_.find(prop_name) != properties_.end();
    }

    std::vector<std::string> Material::getPropertyNames() const {
        std::vector<std::string> names;
        names.reserve(properties_.size());
        for (const auto& pair : properties_) {
            names.push_back(pair.first);
        }
        return names;
    }

} // namespace Core