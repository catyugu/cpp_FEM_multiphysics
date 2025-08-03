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

        // 获取单元几何信息
        const ElementGeometry& getGeometry();
        void update_geometry();
        void setMaterialID(int id) { material_id_ = id; }
        int getMaterialID() const { return material_id_; }
        virtual std::unique_ptr<FEValues> createFEValues(int quad_order) = 0;

    protected:
        int id_;
        std::vector<Node*> nodes_;
        int order_ = 1;
        int material_id_ = 0;
        std::unique_ptr<ElementGeometry> geometry_;

        void set_nodes_internal(const std::vector<Node*>& new_nodes);
    };

} // namespace Core
#endif // ELEMENT_HPP