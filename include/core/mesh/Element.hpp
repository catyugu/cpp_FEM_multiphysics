#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "Node.hpp"
#include <vector>
#include <memory>

namespace Core {

    // Abstract base class for all element types
    class Element {
    public:
        Element(int id) : id_(id) {}
        virtual ~Element() = default;

        // Get the element ID
        int getId() const;

        // Get the nodes connected to this element
        const std::vector<Node*>& getNodes() const;

        // Add a node to the element
        void addNode(Node* node);

        // Pure virtual function to get the number of nodes
        virtual size_t getNumNodes() const = 0;

        // Pure virtual function to get the element type name
        virtual const char* getTypeName() const = 0;


    protected:
        int id_;
        std::vector<Node*> nodes_; // Pointers to the nodes of the element
    };




} // namespace Core

#endif // ELEMENT_HPP
