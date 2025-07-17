#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "Node.hpp"
#include <vector>
#include <memory>
#include "core/mesh/ElementGeometry.hpp" // Include the new header

namespace Core {

    // Forward declaration to break circular dependency
    class FEValues;

    class Element {
    public:
        Element(int id);
        virtual ~Element() = default;

        int getId() const;
        const std::vector<Node*>& getNodes() const;
        void addNode(Node* node);

        virtual size_t getNumNodes() const = 0;
        virtual const char* getTypeName() const = 0;
        virtual int getDimension() const = 0; // <-- ADDED: Pure virtual function

        int getOrder() const { return order_; }
        void setOrder(int order) { order_ = order; }

        // The Element's primary new role: creating FEValues calculators.
        std::unique_ptr<FEValues> create_fe_values(int quad_order);
        // Helper to initialize geometry once nodes are added.
        void update_geometry();
    protected:
        int id_;
        std::vector<Node*> nodes_;
        int order_ = 1;

        // Each element now holds its own geometry.
        std::unique_ptr<ElementGeometry> geometry_;


    };

} // namespace Core
#endif // ELEMENT_HPP