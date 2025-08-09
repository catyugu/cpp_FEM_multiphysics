#include "utils/InterpolationUtilities.hpp"
#include "core/mesh/Element.hpp"
#include "core/material/Material.hpp"
#include "physics/PhysicsField.hpp"
#include "core/DOFManager.hpp"

namespace Utils {

std::map<std::string, double> InterpolationUtilities::interpolateAtQuadraturePoint(
    const Core::Element* element,
    const Eigen::VectorXd& shape_values,
    const std::vector<std::string>& variable_names,
    const std::map<std::string, const Physics::PhysicsField*>& physics_fields) {

    std::map<std::string, double> interpolated_values;

    for (const std::string& var_name : variable_names) {
        auto field_it = physics_fields.find(var_name);
        if (field_it == physics_fields.end()) {
            // 如果找不到对应的物理场，使用元素存储的变量值作为常数
            interpolated_values[var_name] = element->getVariableValue(var_name);
            continue;
        }

        const Physics::PhysicsField* field = field_it->second;
        const auto& solution = field->getSolution();
        const auto& nodes = element->getNodes();

        double interpolated_value = 0.0;

        // 使用形函数进行插值：u(ξ) = Σ N_i(ξ) * u_i
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (i < static_cast<size_t>(shape_values.size())) {
                int dof_idx = field->getDofManager()->getEquationIndex(nodes[i]->getId(), field->getVariableName());
                if (dof_idx != -1 && dof_idx < solution.size()) {
                    interpolated_value += shape_values(i) * solution(dof_idx);
                }
            }
        }

        interpolated_values[var_name] = interpolated_value;
    }

    return interpolated_values;
}

double InterpolationUtilities::evaluateMaterialPropertyAtQuadraturePoint(
    const Core::Material& material,
    const std::string& property_name,
    const std::map<std::string, double>& interpolated_variables) {

    // 首先尝试使用新的MaterialProperty系统
    if (material.hasProperty(property_name)) {
        // 使用插值得到的变量值来评估材料属性
        return material.getProperty(property_name, interpolated_variables);
    }

    // 如果没有找到属性，抛出异常
    throw std::runtime_error("Material property '" + property_name + "' not found");
}

} // namespace Utils
