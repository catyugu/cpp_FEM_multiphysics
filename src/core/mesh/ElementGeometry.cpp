#include "core/mesh/ElementGeometry.hpp"
#include "utils/Exceptions.hpp"

namespace Core {

    // Constructor implementation is updated.
    ElementGeometry::ElementGeometry(const std::vector<Node*>& vertex_nodes, int dimension)
        : vertex_nodes_(vertex_nodes), dimension_(dimension) {
        if (vertex_nodes.empty()) {
            throw Exception::ConfigurationException("ElementGeometry cannot be created with zero nodes.");
        }
        // Resize the coordinate matrix based on the *correct* dimension.
        vertex_coords_.resize(dimension_, vertex_nodes.size());

        for (size_t i = 0; i < vertex_nodes.size(); ++i) {
            const auto& coords = vertex_nodes[i]->getCoords();
            if (coords.size() < dimension_) {
                throw Exception::ConfigurationException("Node has fewer coordinates than the element's dimension.");
            }
            // Only copy the relevant coordinates (e.g., x,y for 2D; x for 1D).
            for(int j = 0; j < dimension_; ++j) {
                vertex_coords_(j, i) = coords[j];
            }
        }
    }

} // namespace Core