#include <core/mesh/Element.hpp>
#include <stdexcept>
#include <cmath>

namespace Core {

    // --- Element Base Class ---
    int Element::getId() const {
        return id_;
    }

    const std::vector<Node*>& Element::getNodes() const {
        return nodes_;
    }

    void Element::addNode(Node* node) {
        nodes_.push_back(node);
    }



} // namespace Core
