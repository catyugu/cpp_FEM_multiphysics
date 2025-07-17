#include <core/mesh/Element.hpp>
#include <core/FEValues.hpp>
#include <stdexcept>

namespace Core {

    int Element::getId() const {
        return id_;
    }

    const std::vector<Node*>& Element::getNodes() const {
        return nodes_;
    }

    void Element::addNode(Node* node) {
        nodes_.push_back(node);
    }

    // FIX: Implementation of the const (read-only) version
    const ElementGeometry& Element::getGeometry() const {
        if (!geometry_) {
            throw std::runtime_error("Element geometry was not initialized before use in element ID: " + std::to_string(id_));
        }
        return *geometry_;
    }

    // FIX: Implementation of the non-const (initializing) version
    ElementGeometry& Element::getGeometry() {
        if (!geometry_) {
            geometry_ = std::make_unique<ElementGeometry>(nodes_);
        }
        return *geometry_;
    }

    std::unique_ptr<FEValues> Element::create_fe_values(int order, int quad_order) {
        return std::make_unique<FEValues>(*this, order, quad_order);
    }

} // namespace Core