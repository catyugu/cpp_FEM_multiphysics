// include/utils/NodeReordering.hpp

#ifndef NODEREORDERING_HPP
#define NODEREORDERING_HPP

#include "core/mesh/Mesh.hpp"
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include "utils/SimpleLogger.hpp"

namespace Utils {

    // A function to compute the RCM permutation for a given mesh
    inline std::vector<int> compute_rcm_permutation(const Core::Mesh& mesh) {
        auto& logger = Logger::instance();
        logger.info("Computing Reverse Cuthill-McKee (RCM) node permutation...");

        // Define the graph type using Boost Graph Library
        using Graph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;
        const auto& nodes = mesh.getNodes();
        int num_nodes = nodes.size();

        // Create the graph from the mesh connectivity
        Graph g(num_nodes);
        for (const auto& elem : mesh.getElements()) {
            const auto& elem_nodes = elem->getNodes();
            for (size_t i = 0; i < elem_nodes.size(); ++i) {
                for (size_t j = i + 1; j < elem_nodes.size(); ++j) {
                    boost::add_edge(elem_nodes[i]->getId(), elem_nodes[j]->getId(), g);
                }
            }
        }

        // The permutation vector from Boost
        std::vector<int> perm(num_nodes);

        // Call the RCM algorithm
        boost::cuthill_mckee_ordering(g, perm.rbegin()); // .rbegin() makes it Reverse Cuthill-McKee

        logger.info("RCM permutation computed.");

        // We need to return an inverse permutation map for easy lookup: new_id -> old_id
        std::vector<int> inv_perm(num_nodes);
        for(int i = 0; i < num_nodes; ++i) {
            inv_perm[perm[i]] = i;
        }

        return inv_perm; // Returns a map where index = new_node_id, value = old_node_id
    }

} // namespace Utils

#endif // NODEREORDERING_HPP