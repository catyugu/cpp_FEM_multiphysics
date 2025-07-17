#include <core/mesh/Element.hpp>
#include <core/FEValues.hpp> // Include the new header
#include <stdexcept>

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

    void Element::update_geometry() {
        if (!nodes_.empty()) {
            // Pass the element's intrinsic dimension to the geometry constructor.
            geometry_ = std::make_unique<ElementGeometry>(nodes_, this->getDimension());
        }
    }

    std::unique_ptr<FEValues> Element::create_fe_values(int quad_order) {
        if (!geometry_) {
            // Lazy initialization of geometry if not already done.
            update_geometry();
        }
        // The element order for the FEValues is taken from the element itself.
        return std::make_unique<FEValues>(*geometry_, this->order_, quad_order);
    }

} // namespace Core