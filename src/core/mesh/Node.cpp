#include <core/mesh/Node.hpp>

namespace Core {

    Node::Node(int id, double x, double y, double z) : id_(id) {
        coords_.push_back(x);
        coords_.push_back(y);
        coords_.push_back(z);
    }

    int Node::getId() const {
        return id_;
    }

    const std::vector<double>& Node::getCoords() const {
        return coords_;
    }

} // namespace Core
