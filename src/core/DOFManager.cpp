#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Node.hpp" // Required to get node details from the mesh

namespace Core {

    DOFManager::DOFManager(Mesh& mesh) : mesh_(mesh), num_equations_(0) {}

    void DOFManager::registerVariable(const std::string& var_name) {
        variable_names_.push_back(var_name);
        SimpleLogger::Logger::instance().info("DOFManager: Registered variable '", var_name, "'.");
    }

    void DOFManager::build() {
        auto& logger = SimpleLogger::Logger::instance();
        logger.info("DOFManager: Building DOF map...");
        dof_map_.clear();
        int equation_counter = 0;

        // A simple DOF mapping: iterate through nodes, then variables
        for (const auto& node : mesh_.getNodes()) {
            for (size_t i = 0; i < variable_names_.size(); ++i) {
                dof_map_[{node->getId(), static_cast<int>(i)}] = equation_counter++;
            }
        }
        num_equations_ = equation_counter;
        logger.info("DOFManager: Built map with ", num_equations_, " equations.");
    }

    int DOFManager::getEquationIndex(int node_id, const std::string& var_name) const {
        int var_idx = -1;
        for (size_t i = 0; i < variable_names_.size(); ++i) {
            if (variable_names_[i] == var_name) {
                var_idx = static_cast<int>(i);
                break;
            }
        }
        if (var_idx == -1) {
            SimpleLogger::Logger::instance().error("DOFManager: Variable '", var_name, "' not registered.");
            return -1;
        }

        auto it = dof_map_.find({node_id, var_idx});
        if (it == dof_map_.end()) {
            SimpleLogger::Logger::instance().error("DOFManager: DOF for node ", node_id, " and var '", var_name, "' not found.");
            return -1;
        }
        return it->second;
    }

    size_t DOFManager::getNumEquations() const {
        return num_equations_;
    }

    const std::vector<std::string>& DOFManager::getVariableNames() const {
        return variable_names_;
    }

} // namespace Core
