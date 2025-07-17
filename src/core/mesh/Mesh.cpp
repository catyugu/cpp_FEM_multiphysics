#include <core/mesh/Mesh.hpp>
#include <core/mesh/Element.hpp>
#include <core/mesh/LineElement.hpp>
#include <core/mesh/TriElement.hpp>
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"

namespace Core {

Mesh::~Mesh() {
    for (auto node : nodes_) {
        delete node;
    }
    for (auto element : elements_) {
        delete element;
    }
}

void Mesh::addNode(Node* node) {
    nodes_.push_back(node);
    node_map_[node->getId()] = node;
}

void Mesh::addElement(Element* element) {
    elements_.push_back(element);
    element_map_[element->getId()] = element;
}

Node *Mesh::getNode(int id) const {
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
    auto& logger = Utils::Logger::instance();
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
        elem->update_geometry();
        mesh->addElement(elem);
    }
    logger.info("Created ", mesh->getNodes().size(), " nodes and ", mesh->getElements().size(), " elements.");
    return mesh;
}

Mesh* Mesh::create_uniform_2d_mesh(double width, double height, int nx, int ny) {
    auto& logger = Utils::Logger::instance();
    logger.info("Creating uniform 2D mesh...");
    Mesh* mesh = new Mesh();
    double dx = width / nx;
    double dy = height / ny;
    int node_id = 0;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            mesh->addNode(new Node(node_id++, i * dx, j * dy));
        }
    }
    int elem_id = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int n0 = j * (nx + 1) + i;
            int n1 = n0 + 1;
            int n2 = (j + 1) * (nx + 1) + i;
            int n3 = n2 + 1;
            TriElement* elem1 = new TriElement(elem_id++);
            elem1->addNode(mesh->getNode(n0));
            elem1->addNode(mesh->getNode(n1));
            elem1->addNode(mesh->getNode(n2));
            elem1->update_geometry();
            mesh->addElement(elem1);
            TriElement* elem2 = new TriElement(elem_id++);
            elem2->addNode(mesh->getNode(n1));
            elem2->addNode(mesh->getNode(n3));
            elem2->addNode(mesh->getNode(n2));
            elem2->update_geometry();
            mesh->addElement(elem2);
        }
    }
    logger.info("Created ", mesh->getNodes().size(), " nodes and ", mesh->getElements().size(), " elements.");
    return mesh;
}

Mesh* Mesh::create_uniform_3d_mesh(double width, double height, double depth, int nx, int ny, int nz) {
    auto& logger = Utils::Logger::instance();
    logger.info("Creating uniform 3D mesh...");
    Mesh* mesh = new Mesh();
    double dx = width / nx;
    double dy = height / ny;
    double dz = depth / nz;

    // Create nodes
    int node_id = 0;
    for (int k = 0; k <= nz; ++k) {
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                mesh->addNode(new Node(node_id++, i * dx, j * dy, k * dz));
            }
        }
    }

    // Create tetrahedral elements by splitting each hexahedral cell into 6 tetrahedrons.
    // The node ordering is crucial here to ensure positive Jacobian determinants.
    int elem_id = 0;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                // Local node indices for the current hexahedral cell
                int n[8];
                n[0] = k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
                n[1] = n[0] + 1;
                n[2] = n[0] + (nx + 1);
                n[3] = n[1] + (nx + 1);
                n[4] = n[0] + (nx + 1) * (ny + 1);
                n[5] = n[1] + (nx + 1) * (ny + 1);
                n[6] = n[2] + (nx + 1) * (ny + 1);
                n[7] = n[3] + (nx + 1) * (ny + 1);

                // Define the 6 tetrahedra with corrected node ordering for positive Jacobian
                auto add_tet = [&](int i0, int i1, int i2, int i3) {
                    auto* tet = new TetElement(elem_id++);
                    tet->addNode(mesh->getNode(n[i0]));
                    tet->addNode(mesh->getNode(n[i1]));
                    tet->addNode(mesh->getNode(n[i2]));
                    tet->addNode(mesh->getNode(n[i3]));
                    tet->update_geometry();
                    mesh->addElement(tet);
                };

                add_tet(0, 1, 3, 7); // Tet 1: Correct
                add_tet(0, 1, 7, 5); // Tet 2: Corrected (swapped 5 and 7 from original)
                add_tet(0, 2, 7, 3); // Tet 3: Corrected (swapped 3 and 7 from original)
                add_tet(0, 2, 6, 7); // Tet 4: Correct
                add_tet(0, 4, 5, 7); // Tet 5: Correct
                add_tet(0, 4, 7, 6); // Tet 6: Corrected (swapped 6 and 7 from original)
            }
        }
    }

    logger.info("Created ", mesh->getNodes().size(), " nodes and ", mesh->getElements().size(), " elements.");
    return mesh;
}


} // namespace Core