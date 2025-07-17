#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "Node.hpp"
#include <vector>
#include <memory>

namespace Core {
    class FEValues;

    class Element {
    public:
        Element(int id) : id_(id) {}
        virtual ~Element() = default;

        int getId() const;
        const std::vector<Node*>& getNodes() const;
        void addNode(Node* node);

        virtual size_t getNumNodes() const = 0;
        virtual const char* getTypeName() const = 0;

        int getOrder() const { return order_; }
        void setOrder(int order) { order_ = order; }

        std::unique_ptr<FEValues> create_fe_values(int order, int quad_order);

    protected:
        int id_;
        std::vector<Node*> nodes_;
        int order_ = 1;
    };

} // namespace Core
#endif // ELEMENT_HPP