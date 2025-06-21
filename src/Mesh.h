#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <Eigen/Dense>

// A simple structure to hold the coordinates of a single node.
struct Node {
    int id;
    Eigen::Vector3d coordinates;
};

// A simple structure to hold the connectivity of a single element.
struct Element {
    int id;
    Eigen::VectorXi node_ids;

    Element(int element_vertice_num, int id) {
        node_ids = Eigen::VectorXi::Zero(element_vertice_num);
        this->id = id;
    }
};

// The main Mesh class.
class Mesh {
public:
    int element_vertice_num = 4;
    std::vector<Node> nodes;
    std::vector<Element> elements;

    // Constructor now defaults to 4 vertices for 3D tetrahedral elements.
    explicit Mesh(int element_vertice_num = 4) {
        this->element_vertice_num = element_vertice_num;
    };

    // Updated public method to route to the correct mesh parser.
    bool load_mesh(const std::string& filename);

private:
    // Helper function for PLY format.
    bool load_mesh_from_ply(const std::string& filename);
    // New helper function for Nastran BDF format.
    bool load_mesh_from_nastran(const std::string& filename);
    // Helper function for TXT format.
    bool load_mesh_from_txt(const std::string& filename);
};

#endif // MESH_H