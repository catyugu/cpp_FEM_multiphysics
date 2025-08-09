#include "core/material/Variable.hpp"
#include <stdexcept>

namespace Core {

Variable::Variable(const std::string& name, double default_value, const std::string& description)
    : name_(name), default_value_(default_value), description_(description) {
    if (name.empty()) {
        throw std::invalid_argument("Variable name cannot be empty");
    }
}

} // namespace Core
