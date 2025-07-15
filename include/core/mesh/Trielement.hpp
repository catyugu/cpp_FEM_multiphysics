#ifndef TRIELEMENT_HPP
#define TRIELEMENT_HPP

#include "Element.hpp"
#include <Eigen/Dense>

namespace Core {

    /**
     * @class TriElement
     * @brief Represents a 2D, 3-node, linear triangular element.
     */
    class TriElement : public Element {
    public:
        explicit TriElement(int id);

        size_t getNumNodes() const override;
        const char* getTypeName() const override;

        // Calculates the area of the triangular element using the shoelace formula.
        double getArea() const;

        // Calculates the B-matrix (strain-displacement matrix) for this element.
        // For 2D heat transfer, this matrix relates the temperature gradients to the nodal temperatures.
        // [ T,x ] = [b1, b2, b3] * { T1 }
        // [ T,y ]   [c1, c2, c3]   { T2 }
        //                          { T3 }
        // where B = (1/2A) * [[y2-y3, y3-y1, y1-y2], [x3-x2, x1-x3, x2-x1]]
        Eigen::Matrix<double, 2, 3> getBMatrix() const;
    };

} // namespace Core

#endif // TRIELEMENT_HPP