#include "core/material/MaterialProperty.hpp"
#include "core/mesh/Element.hpp"
#include <stdexcept>

namespace Core {

MaterialProperty::MaterialProperty(const std::string& name, double value)
    : name_(name), is_constant_(true), constant_value_(value) {
    if (name.empty()) {
        throw std::invalid_argument("Property name cannot be empty");
    }
}

MaterialProperty::MaterialProperty(const std::string& name,
                                 PropertyFunction function,
                                 const std::vector<std::string>& dependent_variables)
    : name_(name), is_constant_(false), constant_value_(0.0),
      function_(std::move(function)), dependent_variables_(dependent_variables) {
    if (name.empty()) {
        throw std::invalid_argument("Property name cannot be empty");
    }
    if (!function_) {
        throw std::invalid_argument("Property function cannot be null");
    }
}

double MaterialProperty::evaluate(const Element* element) const {
    if (is_constant_) {
        return constant_value_;
    }

    if (!element) {
        throw std::invalid_argument("Element cannot be null for variable-dependent property");
    }

    // 从元素获取所有依赖变量的值
    std::map<std::string, double> variable_values;
    for (const std::string& var_name : dependent_variables_) {
        variable_values[var_name] = element->getVariableValue(var_name);
    }

    return function_(variable_values);
}

double MaterialProperty::evaluate(const std::map<std::string, double>& variable_values) const {
    if (is_constant_) {
        return constant_value_;
    }

    return function_(variable_values);
}

void MaterialProperty::setConstantValue(double value) {
    is_constant_ = true;
    constant_value_ = value;
    function_ = nullptr;
    dependent_variables_.clear();
}

void MaterialProperty::setFunction(PropertyFunction function,
                                 const std::vector<std::string>& dependent_variables) {
    if (!function) {
        throw std::invalid_argument("Property function cannot be null");
    }

    is_constant_ = false;
    constant_value_ = 0.0;
    function_ = std::move(function);
    dependent_variables_ = dependent_variables;
}

} // namespace Core
