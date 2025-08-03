#ifndef TETELEMENT_HPP
#define TETELEMENT_HPP

#include <core/FEValues.hpp>
#include <core/ReferenceElement.hpp>

#include "Element.hpp"
#include <Eigen/Dense>

namespace Core {

    /**
     * @class TetElement
     * @brief Represents a 3D, 4-node, linear tetrahedral element.
     */
    class TetElement : public Element {
    public:
        explicit TetElement(int id);

        size_t getNumNodes() const override;
        const char* getTypeName() const override;

        int getDimension() const override { return 3;}

        // Calculates the volume of the tetrahedral element.
        double getVolume() const;

        // Calculates the B-matrix (strain-displacement matrix) for this element.
        // For 3D heat transfer, this 3x4 matrix relates the temperature gradients
        // (in x, y, z) to the 4 nodal temperatures.
        Eigen::Matrix<double, 3, 4> getBMatrix() const;

        std::unique_ptr<FEValues> createFEValues(int quad_order) override;
    };

} // namespace Core

#endif // TETELEMENT_HPP
