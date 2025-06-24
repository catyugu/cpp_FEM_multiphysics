#include "core/Mesh.hpp"
#include "utils/SimpleLogger.hpp"

namespace Core {

    Mesh::~Mesh() {
        for (auto node : nodes_) {
            delete node;
        }
        for (auto element : elements_) {
            delete element;
        }
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

    const std::vector<Element*>& Mesh::getElements() const {
        return elements_;
    }

    const std::vector<Node*>& Mesh::getNodes() const {
        return nodes_;
    }


    Mesh* Mesh::create_uniform_1d_mesh(double length, int num_elements) {
        auto& logger = SimpleLogger::Logger::instance();
        logger.info("Creating uniform 1D mesh...");
        logger.info("Length: ", length, ", Elements: ", num_elements);

        Mesh* mesh = new Mesh();
        double h = length / num_elements;

        // Create nodes
        for (int i = 0; i <= num_elements; ++i) {
            mesh->addNode(new Node(i, i * h));
        }

        // Create elements
        for (int i = 0; i < num_elements; ++i) {
            Element* elem = new LineElement(i);
            elem->addNode(mesh->getNode(i));
            elem->addNode(mesh->getNode(i + 1));
            mesh->addElement(elem);
        }

        logger.info("Created ", mesh->getNodes().size(), " nodes and ", mesh->getElements().size(), " elements.");
        return mesh;
    }


} // namespace Core
