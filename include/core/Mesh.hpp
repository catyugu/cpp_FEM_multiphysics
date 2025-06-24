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

        // Add a node to the mesh
        void addNode(Node* node);

        // Add an element to the mesh
        void addElement(Element* element);

        // Get a node by its ID
        Node* getNode(int id) const;

        // Get an element by its ID
        Element* getElement(int id) const;

        // Get all elements
        const std::vector<Element*>& getElements() const;

        // Get all nodes
        const std::vector<Node*>& getNodes() const;


        // A simple method to create a uniform 1D mesh
        static Mesh* create_uniform_1d_mesh(double length, int num_elements);


    private:
        std::vector<Node*> nodes_;
        std::vector<Element*> elements_;

        // Maps for quick access
        std::map<int, Node*> node_map_;
        std::map<int, Element*> element_map_;
    };

} // namespace Core

#endif // MESH_HPP
