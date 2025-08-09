#include "core/material/VariableManager.hpp"
#include <stdexcept>

namespace Core {

VariableManager& VariableManager::getInstance() {
    static VariableManager instance;
    return instance;
}

const Variable& VariableManager::registerVariable(const std::string& name,
                                                 double default_value,
                                                 const std::string& description) {
    if (name.empty()) {
        throw std::invalid_argument("Variable name cannot be empty");
    }

    if (variables_.find(name) != variables_.end()) {
        throw std::invalid_argument("Variable '" + name + "' already exists");
    }

    auto variable = std::make_unique<Variable>(name, default_value, description);
    const Variable& ref = *variable;
    variables_[name] = std::move(variable);

    return ref;
}

const Variable& VariableManager::getVariable(const std::string& name) const {
    auto it = variables_.find(name);
    if (it == variables_.end()) {
        throw std::out_of_range("Variable '" + name + "' not found");
    }
    return *it->second;
}

bool VariableManager::hasVariable(const std::string& name) const {
    return variables_.find(name) != variables_.end();
}

std::vector<std::string> VariableManager::getAllVariableNames() const {
    std::vector<std::string> names;
    names.reserve(variables_.size());

    for (const auto& pair : variables_) {
        names.push_back(pair.first);
    }

    return names;
}

void VariableManager::clear() {
    variables_.clear();
}

} // namespace Core
