#include <core/mesh/TetElement.hpp>
#include <stdexcept>
#include "utils/ShapeFunctions.hpp"

namespace Core {

TetElement::TetElement(int id) : Element(id) {}

size_t TetElement::getNumNodes() const {
    if (order_ > 0 && order_ <= Utils::MAX_TET_ORDER_SUPPORTED) {
        return (order_ + 1) * (order_ + 2) * (order_ + 3) / 6;
    }
    throw std::runtime_error("Unsupported order " + std::to_string(order_) + " for TetElement. Maximum supported order is " + std::to_string(Utils::MAX_TET_ORDER_SUPPORTED) + ".");
}

const char* TetElement::getTypeName() const {
    return "TetElement";
}

double TetElement::getVolume() const {
    if (nodes_.size() != 4) {
        throw std::runtime_error("TetElement must have exactly 4 nodes to calculate volume.");
    }
    const auto& p1 = nodes_[0]->getCoords();
    const auto& p2 = nodes_[1]->getCoords();
    const auto& p3 = nodes_[2]->getCoords();
    const auto& p4 = nodes_[3]->getCoords();

    Eigen::Matrix4d V_mat;
    V_mat << 1, p1[0], p1[1], p1[2],
             1, p2[0], p2[1], p2[2],
             1, p3[0], p3[1], p3[2],
             1, p4[0], p4[1], p4[2];

    // REMOVED: std::abs() to return the signed volume. This is crucial for detecting inverted elements.
    return V_mat.determinant() / 6.0;
}

Eigen::Matrix<double, 3, 4> TetElement::getBMatrix() const {
    if (nodes_.size() != 4) {
        throw std::runtime_error("TetElement must have exactly 4 nodes to calculate B-matrix.");
    }
    const auto& p1 = nodes_[0]->getCoords();
    const auto& p2 = nodes_[1]->getCoords();
    const auto& p3 = nodes_[2]->getCoords();
    const auto& p4 = nodes_[3]->getCoords();

    double x1 = p1[0], y1 = p1[1], z1 = p1[2];
    double x2 = p2[0], y2 = p2[1], z2 = p2[2];
    double x3 = p3[0], y3 = p3[1], z3 = p3[2];
    double x4 = p4[0], y4 = p4[1], z4 = p4[2];
    Eigen::Matrix4d C;
    C << 1, x1, y1, z1,
         1, x2, y2, z2,
         1, x3, y3, z3,
         1, x4, y4, z4;

    // CHECK: Use std::abs() here for the tolerance check, to detect truly degenerate (zero volume) elements.
    if (std::abs(C.determinant()) < 1e-12) {
        throw std::runtime_error("Element " + std::to_string(id_) + " has zero or near-zero volume, cannot calculate B-matrix.");
    }

    Eigen::Matrix4d C_inv = C.inverse();

    Eigen::Matrix<double, 3, 4> B;
    B.row(0) = C_inv.row(1);
    B.row(1) = C_inv.row(2);
    B.row(2) = C_inv.row(3);

    return B;
}

} // namespace Core