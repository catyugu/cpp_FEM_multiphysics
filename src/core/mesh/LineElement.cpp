#include <core/mesh/LineElement.hpp>
#include <stdexcept>
#include <cmath>
#include "utils/ShapeFunctions.hpp"


namespace Core {

    // --- LineElement ---
    LineElement::LineElement(int id) : Element(id) {}

    size_t LineElement::getNumNodes() const {
        if (order_ > 0 && order_ <= Utils::MAX_LINE_ORDER_SUPPORTED) {
            // The number of nodes for a line element of order p is p+1
            return order_ + 1;
        }
        throw std::runtime_error("Unsupported order " + std::to_string(order_) + " for LineElement. Maximum supported order is " + std::to_string(Utils::MAX_LINE_ORDER_SUPPORTED) + ".");
    }

    const char* LineElement::getTypeName() const {
        return "LineElement";
    }

    double LineElement::getLength() const {
        if (nodes_.size() != 2) {
            throw std::runtime_error("LineElement must have exactly 2 nodes to calculate length.");
        }
        const auto& p1 = nodes_[0]->getCoords();
        const auto& p2 = nodes_[1]->getCoords();
        return std::sqrt(std::pow(p2[0] - p1[0], 2) +
                         std::pow(p2[1] - p1[1], 2) +
                         std::pow(p2[2] - p1[2], 2));
    }
    std::unique_ptr<FEValues> LineElement::createFEValues(int quad_order) {
        const auto& ref_data = ReferenceElementCache::get(getTypeName(), getNodes().size(), getOrder(), quad_order);
        return std::make_unique<FEValues>(getGeometry(), getOrder(), ref_data);
    }

}