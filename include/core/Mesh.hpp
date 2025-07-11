#ifndef MESH_HPP
#define MESH_HPP

#include "Node.hpp"
#include "Element.hpp"
#include <vector>
#include <string>
#include <map>

namespace Core {

    class Mesh {
    public:
        ~Mesh();

        void addNode(Node* node);
        void addElement(Element* element);

        Node *getNode(int id) const;
        Element* getElement(int id) const;

        const std::vector<Element*>& getElements() const;
        const std::vector<Node*>& getNodes() const;

        // --- Mesh Generators ---
        static Mesh* create_uniform_1d_mesh(double length, int num_elements);
        static Mesh* create_uniform_2d_mesh(double width, double height, int nx, int ny);
        static Mesh* create_uniform_3d_mesh(double width, double height, double depth, int nx, int ny, int nz);

    private:
        std::vector<Node*> nodes_;
        std::vector<Element*> elements_;

        // Maps for quick access
        std::map<int, Node*> node_map_;
        std::map<int, Element*> element_map_;
    };

} // namespace Core

#endif // MESH_HPP
