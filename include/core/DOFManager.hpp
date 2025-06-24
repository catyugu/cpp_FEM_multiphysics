#ifndef DOFMANAGER_HPP
#define DOFMANAGER_HPP

#include "Mesh.hpp"
#include <vector>
#include <map>
#include <string>

namespace Core {

    class DOFManager {
    public:
        // Constructor
        DOFManager(Mesh& mesh);

        // Register a degree of freedom for all nodes (e.g., "Temperature")
        void registerVariable(const std::string& var_name);

        // Build the mapping from (node_id, var_name) to a global equation index
        void build();

        // Get the equation index for a given node and variable
        int getEquationIndex(int node_id, const std::string& var_name) const;

        // Get the total number of equations (DOFs)
        size_t getNumEquations() const;

        // Get a list of all registered variable names
        const std::vector<std::string>& getVariableNames() const;

    private:
        Mesh& mesh_;
        std::vector<std::string> variable_names_;
        // Maps a node ID and variable index to a global equation number
        std::map<std::pair<int, int>, int> dof_map_;
        size_t num_equations_;
    };

} // namespace Core

#endif // DOFMANAGER_HPP
