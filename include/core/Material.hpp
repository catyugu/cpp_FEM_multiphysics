#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <string>
#include <map>
#include <any>
#include <functional> // Required for std::function
#include "utils/SimpleLogger.hpp"

namespace Core {

    class Material {
    public:
        Material(int id, const std::string& name);
        using PropertyFunction = std::function<double(const std::map<std::string, double>&)>;
        void setProperty(const std::string& prop_name, double value);
        void setProperty(const std::string& prop_name, PropertyFunction func);
        double getProperty(const std::string& prop_name, const std::map<std::string, double>& field_values = {}) const;

        const std::string& getName() const;
        int getID() const { return id_; }
    private:
        std::string name_;
        std::map<std::string, std::any> properties_;
        int id_;
    };

} // namespace Core

#endif // MATERIAL_HPP