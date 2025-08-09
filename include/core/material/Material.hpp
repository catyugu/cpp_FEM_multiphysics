#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <string>
#include <map>
#include <any>
#include <functional>
#include <vector>
#include "utils/SimpleLogger.hpp"

// Forward declarations
namespace Core {
    class MaterialProperty;
    class Element;
}

namespace Core {

    class Material {
    public:
        Material(int id, std::string  name);
        using PropertyFunction = std::function<double(const std::map<std::string, double>&)>;

        void setProperty(const std::string& prop_name, double value);
        void setProperty(const std::string& prop_name, PropertyFunction func);
        double getProperty(const std::string& prop_name, const std::map<std::string, double>& field_values = {}) const;

        const std::string& getName() const;
        int getID() const { return id_; }

        // ========= 新增：扩展的材料属性管理功能 =========

        /**
         * @brief Set a MaterialProperty object for this material
         * @param prop_name Name of the property
         * @param material_prop MaterialProperty object containing the property definition
         */
        void setMaterialProperty(const std::string& prop_name,
                               const MaterialProperty& material_prop);

        /**
         * @brief Get property value by evaluating it with current element variable values
         * @param prop_name Name of the property to retrieve
         * @param element Element containing current variable values
         * @return Evaluated property value
         */
        double getPropertyFromElement(const std::string& prop_name,
                                    const Element* element) const;

        /**
         * @brief Get property value by evaluating it with interpolated values at a quadrature point
         * @param prop_name Name of the property to retrieve
         * @param interpolated_values Map of variable name to interpolated value at the quadrature point
         * @return Evaluated property value at the quadrature point
         */
        double getPropertyAtQuadraturePoint(const std::string& prop_name,
                                           const std::map<std::string, double>& interpolated_values) const;

        /**
         * @brief Check if a property exists
         * @param prop_name Name of the property to check
         * @return True if property exists, false otherwise
         */
        bool hasProperty(const std::string& prop_name) const;

        /**
         * @brief Get all property names
         * @return Vector of property names
         */
        std::vector<std::string> getPropertyNames() const;

    private:
        std::string name_;
        std::map<std::string, std::any> properties_;
        int id_;
    };

} // namespace Core

#endif // MATERIAL_HPP