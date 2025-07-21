#include "core/DOFManager.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/Element.hpp>

namespace Core {

DOFManager::DOFManager(Mesh& mesh) : mesh_(mesh), num_equations_(0) {}

void DOFManager::registerVariable(const std::string& var_name, int num_components) {
    variable_names_.push_back(var_name);
    variable_components_[var_name] = num_components;
}
void DOFManager::build(const std::map<std::string, int>& field_orders) {
    auto& logger = Utils::Logger::instance();
    logger.info("DOFManager: Building advanced DOF map...");

    vertex_dof_map_.clear();
    edge_dof_map_.clear();
    int equation_counter = 0;

    // --- Pass 1: Assign DOFs to all vertex nodes ---
    for (const auto& node : mesh_.getNodes()) {
        for (size_t i = 0; i < variable_names_.size(); ++i){
            int num_components = variable_components_.at(variable_names_[i]);
            for (int c = 0; c < num_components; ++c) {
                vertex_dof_map_[{node->getId(), i}] = equation_counter;
                equation_counter++;
            }
        }
    }

    // --- Pass 2: Assign DOFs to higher-order nodes (edges) if needed ---
    for (const auto& elem_ptr : mesh_.getElements()) {
        for (size_t var_idx = 0; var_idx < variable_names_.size(); ++var_idx) {
            const auto& var_name = variable_names_[var_idx];
            if (field_orders.count(var_name) && field_orders.at(var_name) > 1) {
                // This logic is for quadratic elements (order 2), which have nodes on edge midpoints.
                // It can be extended for cubic+ elements (which also have face/volume nodes).
                const auto& nodes = elem_ptr->getNodes();
                for (size_t i = 0; i < nodes.size(); ++i) {
                    for (size_t j = i + 1; j < nodes.size(); ++j) {
                        std::vector<int> edge_nodes = {nodes[i]->getId(), nodes[j]->getId()};
                        // Sorting guarantees a unique key for each edge regardless of node order
                        std::sort(edge_nodes.begin(), edge_nodes.end());

                        EdgeDofKey key = {edge_nodes, static_cast<int>(var_idx)};
                        if (edge_dof_map_.find(key) == edge_dof_map_.end()) {
                            edge_dof_map_[key] = equation_counter++;
                        }
                    }
                }
            }
        }
    }

    num_equations_ = equation_counter;
    logger.info("DOFManager: Built map with ", num_equations_, " total equations.");
}

int DOFManager::getEquationIndex(int node_id, const std::string& var_name) const {
    int var_idx = -1;
    for (size_t i = 0; i < variable_names_.size(); ++i) if (variable_names_[i] == var_name) var_idx = i;
    if (var_idx == -1) return -1;

    auto it = vertex_dof_map_.find({node_id, var_idx});
    return (it != vertex_dof_map_.end()) ? it->second : -1;
}

int DOFManager::getEdgeEquationIndex(std::vector<int> node_ids, const std::string& var_name) const {
    int var_idx = -1;
    for (size_t i = 0; i < variable_names_.size(); ++i) if (variable_names_[i] == var_name) var_idx = i;
    if (var_idx == -1) return -1;

    std::sort(node_ids.begin(), node_ids.end());
    auto it = edge_dof_map_.find({node_ids, var_idx});
    return (it != edge_dof_map_.end()) ? it->second : -1;
}

size_t DOFManager::getNumEquations() const { return num_equations_; }
const std::vector<std::string>& DOFManager::getVariableNames() const { return variable_names_; }

} // namespace Core