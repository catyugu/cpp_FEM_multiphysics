#ifndef DOFMANAGER_HPP
#define DOFMANAGER_HPP

#include "mesh/Mesh.hpp"
#include <vector>
#include <map>
#include <string>
#include <utility> // Required for std::pair
#include <algorithm> // Required for std::sort

namespace Core {

    // A key for DOFs on vertices (Node ID, Variable Index)
    using VertexDofKey = std::pair<int, int>;

    // A key for higher-order DOFs on edges (Sorted list of Node IDs, Variable Index)
    using EdgeDofKey = std::pair<std::vector<int>, int>;

    class DOFManager {
    public:
        DOFManager(Mesh& mesh);

        void registerVariable(const std::string& var_name);

        // The build method now needs to know the field's requested order
        void build(const std::map<std::string, int>& field_orders);

        // Get equation index for a standard vertex DOF
        int getEquationIndex(int node_id, const std::string& var_name) const;

        // Get equation index for a higher-order edge DOF
        int getEdgeEquationIndex(std::vector<int> node_ids, const std::string& var_name) const;

        size_t getNumEquations() const;
        const std::vector<std::string>& getVariableNames() const;

    private:
        Mesh& mesh_;
        std::vector<std::string> variable_names_;

        // Map for vertex-based DOFs
        std::map<VertexDofKey, int> vertex_dof_map_;

        // New map for edge-based (higher-order) DOFs
        std::map<EdgeDofKey, int> edge_dof_map_;

        size_t num_equations_;
    };

} // namespace Core

#endif // DOFMANAGER_HPP