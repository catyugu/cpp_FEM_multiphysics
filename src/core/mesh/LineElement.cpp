//
// Created by HUAWEI on 2025/7/15.
//

#include <stdexcept>
#include <core/mesh/LineElement.hpp>
#include <cmath>


namespace Core {

    // --- LineElement ---
    LineElement::LineElement(int id) : Element(id) {}

    size_t LineElement::getNumNodes() const {
        return 2;
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

}