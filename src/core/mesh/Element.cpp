#include <core/mesh/Element.hpp>
#include <core/FEValues.hpp>
#include <stdexcept>
#include <algorithm> // Required for std::copy, etc.

namespace Core {

    Element::Element(int id) : id_(id), order_(1), geometry_(nullptr) {}

    int Element::getId() const {
        return id_;
    }

    const std::vector<Node*>& Element::getNodes() const {
        return nodes_;
    }

    void Element::addNode(Node* node) {
        nodes_.push_back(node);
    }

    // NEW: Implementation of set_nodes_internal
    void Element::set_nodes_internal(const std::vector<Node*>& new_nodes) {
        nodes_.clear(); // Clear existing nodes
        nodes_.insert(nodes_.begin(), new_nodes.begin(), new_nodes.end()); // Copy new nodes
    }

    void Element::update_geometry() {
        if (!nodes_.empty()) {
            geometry_ = std::make_unique<ElementGeometry>(nodes_, this->getDimension());
        }
    }

    std::unique_ptr<FEValues> Element::create_fe_values(int quad_order) {
        if (!geometry_) {
            update_geometry();
        }
        return std::make_unique<FEValues>(*geometry_, this->order_, quad_order);
    }

} // namespace Core