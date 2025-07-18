#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "Node.hpp"
#include <vector>
#include <memory>
#include "core/mesh/ElementGeometry.hpp"

// Forward declaration of IO::Importer to declare it as a friend
namespace IO {
    class Importer;
}

namespace Core {

    class FEValues; // Forward declaration

    class Element {
        // NEW: Declare IO::Importer as a friend class
        friend class IO::Importer;

    public:
        Element(int id);
        virtual ~Element() = default;

        int getId() const;
        const std::vector<Node*>& getNodes() const;
        void addNode(Node* node);

        virtual size_t getNumNodes() const = 0;
        virtual const char* getTypeName() const = 0;
        virtual int getDimension() const = 0;

        int getOrder() const { return order_; }
        void setOrder(int order) { order_ = order; }

        std::unique_ptr<FEValues> create_fe_values(int quad_order);
        void update_geometry();
    protected:
        int id_;
        std::vector<Node*> nodes_; // This is the mutable list
        int order_ = 1;
        std::unique_ptr<ElementGeometry> geometry_;

        // Protected method to allow derived classes/friends to set the entire node list
        void set_nodes_internal(const std::vector<Node*>& new_nodes);
    };

} // namespace Core
#endif // ELEMENT_HPP