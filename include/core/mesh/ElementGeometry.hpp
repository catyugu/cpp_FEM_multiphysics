#ifndef ELEMENTGEOMETRY_HPP
#define ELEMENTGEOMETRY_HPP

#include "Node.hpp"
#include <vector>
#include <Eigen/Dense>

namespace Core {

    class ElementGeometry {
    public:
        explicit ElementGeometry(const std::vector<Node*>& vertex_nodes);

        double getVolume() const; // Example for a Tet element
        std::vector<Node *> ElementGeometry::getNodes() const;

    private:
        std::vector<Node *> vertex_nodes_;
        // Caches for geometric properties
        mutable double volume_;
        mutable bool volume_calculated_ = false;
    };

} // namespace Core

#endif // ELEMENTGEOMETRY_HPP