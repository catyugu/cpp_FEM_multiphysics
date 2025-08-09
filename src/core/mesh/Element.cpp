#include <core/mesh/Element.hpp>
#include <core/FEValues.hpp>
#include <core/material/VariableManager.hpp>
#include <stdexcept>
#include <algorithm>

namespace Core {

    Element::Element(int id) : id_(id), order_(1), geometry_(nullptr) {}

    int Element::getId() const {
        return id_;
    }

    const std::vector<Node*>& Element::getNodes() const {
        return nodes_;
    }

    void Element::addNode(Node* node) {
        nodes_.push_back(node);
    }

    void Element::set_nodes_internal(const std::vector<Node*>& new_nodes) {
        nodes_.clear();
        nodes_.insert(nodes_.begin(), new_nodes.begin(), new_nodes.end());
    }

    void Element::update_geometry() {
        if (!nodes_.empty()) {
            geometry_ = std::make_unique<ElementGeometry>(nodes_, this->getDimension());
        }
    }

    // 新增 getGeometry 的实现
    const ElementGeometry& Element::getGeometry() {
        if (!geometry_) {
            update_geometry();
        }
        return *geometry_;
    }

    // 新增：获取与此单元关联的 FEValues 对象，如果尚未创建则创建并初始化
    FEValues* Element::getFEValues(int quad_order, AnalysisType analysis_type) {
        // 检查缓存是否需要更新 (例如，如果积分阶数或分析类型改变)
        if (!fe_values_ || stored_quad_order_ != quad_order || stored_analysis_type_ != analysis_type) {
            fe_values_ = createFEValues(quad_order); // 调用已有的工厂方法
            fe_values_->setAnalysisType(analysis_type); // 设置分析类型并构建B矩阵
            stored_quad_order_ = quad_order;
            stored_analysis_type_ = analysis_type;
        }
        return fe_values_.get();
    }

    // ========= 新增：变量值存储和管理功能实现 =========
    
    void Element::setVariableValue(const std::string& variable_name, double value) {
        variable_values_[variable_name] = value;
    }
    
    double Element::getVariableValue(const std::string& variable_name) const {
        auto it = variable_values_.find(variable_name);
        if (it != variable_values_.end()) {
            return it->second;
        }
        
        // 如果元素没有设置这个变量值，返回变量管理器中的默认值
        try {
            const auto& var_manager = VariableManager::getInstance();
            if (var_manager.hasVariable(variable_name)) {
                return var_manager.getVariable(variable_name).getDefaultValue();
            }
        } catch (const std::exception&) {
            // 变量管理器中也没有这个变量，返回0作为默认值
        }
        
        return 0.0;
    }
    
    bool Element::hasVariableValue(const std::string& variable_name) const {
        return variable_values_.find(variable_name) != variable_values_.end();
    }
    
    const std::map<std::string, double>& Element::getAllVariableValues() const {
        return variable_values_;
    }
    
    void Element::clearVariableValues() {
        variable_values_.clear();
    }
    
    void Element::updateVariableFromSolution(const std::string& variable_name,
                                           const std::vector<double>& solution_values,
                                           const std::vector<int>& dof_indices) {
        if (dof_indices.empty() || solution_values.empty()) {
            return;
        }
        
        // 对于简单的标量场，我们可以计算元素中心的平均值
        double average_value = 0.0;
        size_t valid_dofs = 0;
        
        for (int dof_index : dof_indices) {
            if (dof_index >= 0 && dof_index < static_cast<int>(solution_values.size())) {
                average_value += solution_values[dof_index];
                valid_dofs++;
            }
        }
        
        if (valid_dofs > 0) {
            average_value /= valid_dofs;
            setVariableValue(variable_name, average_value);
        }
    }

} // namespace Core