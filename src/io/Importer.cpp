#include "io/Importer.hpp"
#include <core/mesh/Mesh.hpp>
#include <core/mesh/Node.hpp>
#include <core/mesh/TriElement.hpp> // For 2D meshes
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

namespace IO {

std::unique_ptr<Core::Mesh> Importer::read_comsol_mphtxt(const std::string& filename) {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Importing COMSOL mesh from file: ", filename);

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw Exception::FileIOException("Failed to open mesh file: " + filename);
    }

    auto mesh = std::make_unique<Core::Mesh>();
    std::string line;
    int element_id_counter = 0;

    // --- Two-Pass Parsing Approach for Robustness ---

    // --- Pass 1: Read all Vertices/Nodes ---
    logger.info("Pass 1: Reading mesh vertices...");
    while (std::getline(file, line)) {
        if (line.find("number of mesh vertices") != std::string::npos) {
            std::stringstream ss(line);
            int num_vertices = 0;
            ss >> num_vertices;

            // Find the data header
            while (std::getline(file, line) && line.find("# Mesh vertex coordinates") == std::string::npos);

            // Read the specified number of vertex lines
            for (int i = 0; i < num_vertices && std::getline(file, line); ++i) {
                std::stringstream data_ss(line);
                double x, y;
                if (data_ss >> x >> y) {
                    int node_id = mesh->getNodes().size();
                    mesh->addNode(new Core::Node(node_id, x, y));
                }
            }
            logger.info("Finished reading ", mesh->getNodes().size(), " vertices.");
            break; // Vertices found and read, exit loop.
        }
    }

    // --- Pass 2: Read all Triangle Elements ---
    file.clear(); // Clear any error flags (like EOF)
    file.seekg(0, std::ios::beg); // Rewind file to the beginning
    logger.info("Pass 2: Reading triangle elements...");

    while (std::getline(file, line)) {
        if (line.find("3 tri # type name") != std::string::npos) {
            int num_elements = 0;
            // Find the line that specifies the number of elements for this block
            while (std::getline(file, line) && line.find("number of elements") == std::string::npos);

            std::stringstream ss(line);
            ss >> num_elements;

            // Find the data header
            while (std::getline(file, line) && line.find("# Elements") == std::string::npos);

            // Read the specified number of element lines
            for (int i = 0; i < num_elements && std::getline(file, line); ++i) {
                std::stringstream data_ss(line);
                int n1, n2, n3;
                if (data_ss >> n1 >> n2 >> n3) {
                    Core::Node* node1 = mesh->getNode(n1);
                    Core::Node* node2 = mesh->getNode(n2);
                    Core::Node* node3 = mesh->getNode(n3);

                    if (node1 && node2 && node3) {
                        auto* tri_elem = new Core::TriElement(element_id_counter++);
                        tri_elem->addNode(node1);
                        tri_elem->addNode(node2);
                        tri_elem->addNode(node3);
                        mesh->addElement(tri_elem);
                    } else {
                        logger.error("Invalid node index found in element definition: ", line);
                    }
                }
            }
            logger.info("Finished reading ", mesh->getElements().size(), " elements.");
            break; // Triangles found and read, exit loop.
        }
    }

    logger.info("Import finished. Total: ", mesh->getNodes().size(), " nodes, ", mesh->getElements().size(), " elements.");

    if (mesh->getNodes().empty() || mesh->getElements().empty()) {
        throw Exception::FileIOException("The importer did not read any valid mesh data. Please check the .mphtxt file format and content.");
    }

    return mesh;
}

} // namespace IO