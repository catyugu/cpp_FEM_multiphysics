#include <core/mesh/TriElement.hpp>
#include <stdexcept>
#include "utils/ShapeFunctions.hpp"

namespace Core {

    TriElement::TriElement(int id) : Element(id) {}

    size_t TriElement::getNumNodes() const {
        if (order_ > 0 && order_ <= Utils::MAX_TRI_ORDER_SUPPORTED) {
            // The number of nodes for a triangular element of order p is (p+1)(p+2)/2
            return (order_ + 1) * (order_ + 2) / 2;
        }
        throw std::runtime_error("Unsupported order " + std::to_string(order_) + " for TriElement. Maximum supported order is " + std::to_string(Utils::MAX_TRI_ORDER_SUPPORTED) + ".");
    }

    const char* TriElement::getTypeName() const {
        return "TriElement";
    }

    double TriElement::getArea() const {
        if (nodes_.size() < 3) {
            throw std::runtime_error("TriElement must have at least 3 nodes to calculate area.");
        }
        const auto& p1 = nodes_[0]->getCoords();
        const auto& p2 = nodes_[1]->getCoords();
        const auto& p3 = nodes_[2]->getCoords();
        // Shoelace formula for area
        return 0.5 * std::abs(p1[0]*(p2[1] - p3[1]) + p2[0]*(p3[1] - p1[1]) + p3[0]*(p1[1] - p2[1]));
    }

    Eigen::Matrix<double, 2, 3> TriElement::getBMatrix() const {
        if (nodes_.size() != 3) {
            throw std::runtime_error("TriElement must have exactly 3 nodes to calculate B-matrix.");
        }
        const auto& p1 = nodes_[0]->getCoords();
        const auto& p2 = nodes_[1]->getCoords();
        const auto& p3 = nodes_[2]->getCoords();

        double x1 = p1[0], y1 = p1[1];
        double x2 = p2[0], y2 = p2[1];
        double x3 = p3[0], y3 = p3[1];

        double area = getArea();
        if (area < 1e-12) {
            throw std::runtime_error("Element " + std::to_string(id_) + " has zero or negative area.");
        }

        Eigen::Matrix<double, 2, 3> B;
        B(0, 0) = y2 - y3; B(0, 1) = y3 - y1; B(0, 2) = y1 - y2;
        B(1, 0) = x3 - x2; B(1, 1) = x1 - x3; B(1, 2) = x2 - x1;

        return B / (2.0 * area);
    }

} // namespace Core