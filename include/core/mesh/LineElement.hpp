//
// Created by HUAWEI on 2025/7/15.
//

#ifndef LINEELEMENT_HPP
#define LINEELEMENT_HPP
#include <core/FEValues.hpp>
#include <core/ReferenceElement.hpp>

#include "Element.hpp"

namespace Core{
    // Concrete implementation for a 1D line element with 2 nodes
    class LineElement : public Element {
    public:
        LineElement(int id);
        int getDimension() const override { return 1;}

        size_t getNumNodes() const override;
        const char* getTypeName() const override;
        double getLength() const;
        std::unique_ptr<FEValues> create_fe_values(int quad_order) override;
    };

}




#endif //LINEELEMENT_HPP
