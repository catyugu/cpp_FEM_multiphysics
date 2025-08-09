#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "Node.hpp"
#include <vector>
#include <memory>
#include <map>
#include <string>
#include "core/mesh/ElementGeometry.hpp"
#include "core/AnalysisTypes.hpp" // 添加分析类型头文件

// Forward declaration of IO::Importer to declare it as a friend
namespace IO {
    class Importer;
}

namespace Core {

    class FEValues; // Forward declaration

    class Element {
        friend class IO::Importer;

    public:
        Element(int id);
        virtual ~Element() = default;

        int getId() const;
        const std::vector<Node*>& getNodes() const;
        void addNode(Node* node);

        virtual size_t getNumNodes() const = 0;
        virtual const char* getTypeName() const = 0;
        virtual int getDimension() const = 0;

        int getOrder() const { return order_; }
        void setOrder(int order) { order_ = order; }

        // 获取单元几何信息
        const ElementGeometry& getGeometry();
        void update_geometry();
        void setMaterialID(int id) { material_id_ = id; }
        int getMaterialID() const { return material_id_; }
        virtual std::unique_ptr<FEValues> createFEValues(int quad_order) = 0;

        // 新增：获取与此单元关联的 FEValues 对象，如果尚未创建则创建并初始化
        FEValues* getFEValues(int quad_order, AnalysisType analysis_type);

        // ========= 新增：变量值存储和管理功能 =========

        /**
         * @brief Set the value of a variable for this element
         * @param variable_name Name of the variable
         * @param value Value to set
         */
        void setVariableValue(const std::string& variable_name, double value);

        /**
         * @brief Get the value of a variable for this element
         * @param variable_name Name of the variable
         * @return The variable value, or default value if not set
         */
        double getVariableValue(const std::string& variable_name) const;

        /**
         * @brief Check if a variable value has been set for this element
         * @param variable_name Name of the variable
         * @return True if variable value exists, false otherwise
         */
        bool hasVariableValue(const std::string& variable_name) const;

        /**
         * @brief Get all variable values for this element
         * @return Map of variable name to value
         */
        const std::map<std::string, double>& getAllVariableValues() const;

        /**
         * @brief Clear all variable values
         */
        void clearVariableValues();

        /**
         * @brief Update variable values from a solution vector
         * @param variable_name Name of the variable to update
         * @param solution_values Vector of solution values
         * @param dof_indices DOF indices for this element
         */
        void updateVariableFromSolution(const std::string& variable_name,
                                      const std::vector<double>& solution_values,
                                      const std::vector<int>& dof_indices);

    protected:
        int id_;
        std::vector<Node*> nodes_;
        int order_ = 1;
        int material_id_ = 0;
        std::unique_ptr<ElementGeometry> geometry_;

        // 新增：缓存的 FEValues 对象和相关状态
        std::unique_ptr<FEValues> fe_values_;
        int stored_quad_order_ = -1;  // 用于检查积分阶数是否变化

        // 新增：变量值存储
        mutable std::map<std::string, double> variable_values_;

        void set_nodes_internal(const std::vector<Node*>& new_nodes);
    };

} // namespace Core
#endif // ELEMENT_HPP

