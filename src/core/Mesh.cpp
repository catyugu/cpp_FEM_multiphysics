#include "core/Mesh.hpp"
#include "core/TriElement.hpp" // Include the new 2D element
#include "utils/SimpleLogger.hpp"

namespace Core {

// Destructor and 1D methods remain the same...
Mesh::~Mesh() {
    for (auto node : nodes_) delete node;
    for (auto element : elements_) delete element;
    SimpleLogger::Logger::instance().info("Mesh cleaned up.");
}
void Mesh::addNode(Node* node) {
    nodes_.push_back(node);
    node_map_[node->getId()] = node;
}
void Mesh::addElement(Element* element) {
    elements_.push_back(element);
    element_map_[element->getId()] = element;
}
Node* Mesh::getNode(int id) const {
    auto it = node_map_.find(id);
    return (it != node_map_.end()) ? it->second : nullptr;
}
Element* Mesh::getElement(int id) const {
    auto it = element_map_.find(id);
    return (it != element_map_.end()) ? it->second : nullptr;
}
const std::vector<Element*>& Mesh::getElements() const { return elements_; }
const std::vector<Node*>& Mesh::getNodes() const { return nodes_; }

Mesh* Mesh::create_uniform_1d_mesh(double length, int num_elements) {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Creating uniform 1D mesh...");
    Mesh* mesh = new Mesh();
    double h = length / num_elements;
    for (int i = 0; i <= num_elements; ++i) {
        mesh->addNode(new Node(i, i * h));
    }
    for (int i = 0; i < num_elements; ++i) {
        Element* elem = new LineElement(i);
        elem->addNode(mesh->getNode(i));
        elem->addNode(mesh->getNode(i + 1));
        mesh->addElement(elem);
    }
    logger.info("Created ", mesh->getNodes().size(), " nodes and ", mesh->getElements().size(), " elements.");
    return mesh;
}

// New 2D Mesh Generator
Mesh* Mesh::create_uniform_2d_mesh(double width, double height, int nx, int ny) {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Creating uniform 2D mesh...");
    Mesh* mesh = new Mesh();
    double dx = width / nx;
    double dy = height / ny;

    // Create nodes
    int node_id = 0;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            mesh->addNode(new Node(node_id++, i * dx, j * dy));
        }
    }

    // Create triangular elements
    int elem_id = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i;
            int n1 = n0 + 1;
            int n2 = (j + 1) * (nx + 1) + i;
            int n3 = n2 + 1;

            // Split each quad into two triangles
            TriElement* elem1 = new TriElement(elem_id++);
            elem1->addNode(mesh->getNode(n0));
            elem1->addNode(mesh->getNode(n1));
            elem1->addNode(mesh->getNode(n2));
            mesh->addElement(elem1);

            TriElement* elem2 = new TriElement(elem_id++);
            elem2->addNode(mesh->getNode(n1));
            elem2->addNode(mesh->getNode(n3));
            elem2->addNode(mesh->getNode(n2));
            mesh->addElement(elem2);
        }
    }

    logger.info("Created ", mesh->getNodes().size(), " nodes and ", mesh->getElements().size(), " elements.");
    return mesh;
}

} // namespace Core
