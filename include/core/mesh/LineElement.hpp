//
// Created by HUAWEI on 2025/7/15.
//

#ifndef LINEELEMENT_HPP
#define LINEELEMENT_HPP
#include "Element.hpp"

namespace Core{
    // Concrete implementation for a 1D line element with 2 nodes
    class LineElement : public Element {
    public:
        LineElement(int id);

        size_t getNumNodes() const override;
        const char* getTypeName() const override;
        double getLength() const;
    };

}




#endif //LINEELEMENT_HPP
