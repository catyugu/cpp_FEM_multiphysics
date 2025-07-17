#include "core/mesh/ElementGeometry.hpp"
#include "core/mesh/Node.hpp"
#include <stdexcept>
#include <vector>

namespace Core {

    ElementGeometry::ElementGeometry(const std::vector<Node*>& vertex_nodes) : vertex_nodes_(vertex_nodes), volume_calculated_(false) {}

    // This is an example implementation for a Tetrahedron. You would expand this
    // with `getArea()` for Triangles, `getLength()` for Lines, etc.
    double ElementGeometry::getVolume() const {
        if (volume_calculated_) {
            return volume_;
        }

        if (vertex_nodes_.size() != 4) {
            // Not a tetrahedron, or at least not a simple one.
            // A more robust implementation would handle other element types.
            return 0.0;
        }

        const auto& p1 = vertex_nodes_[0]->getCoords();
        const auto& p2 = vertex_nodes_[1]->getCoords();
        const auto& p3 = vertex_nodes_[2]->getCoords();
        const auto& p4 = vertex_nodes_[3]->getCoords();

        Eigen::Matrix4d V_mat;
        V_mat << 1, p1[0], p1[1], p1[2],
                 1, p2[0], p2[1], p2[2],
                 1, p3[0], p3[1], p3[2],
                 1, p4[0], p4[1], p4[2];

        volume_ = std::abs(V_mat.determinant()) / 6.0;
        volume_calculated_ = true;
        return volume_;
    }

    std::vector<Node *> ElementGeometry::getNodes() const {
        return vertex_nodes_;
    }


} // namespace Core