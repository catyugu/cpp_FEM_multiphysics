#include "mesh.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm> // For std::remove
#include <map>

#include "SimpleLogger.hpp"

// Ensure the global logger is defined
extern SimpleLogger::Logger& logger;

// Helper function to trim whitespace from a string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (std::string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}


// Method to load a mesh file, now routing to PLY or BDF parsers.
bool Mesh::load_mesh(const std::string& filename) {
    logger.info("Loading mesh file: " + filename);
    std::string extension = filename.substr(filename.find_last_of('.') + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    if (extension == "ply") {
        // Attempt to determine if it's PLY or BDF based on content
        std::ifstream file(filename);
        std::string first_line;
        std::getline(file, first_line);
        if (first_line.find("ply") != std::string::npos) {
            logger.info("Detected PLY file format.");
            return load_mesh_from_ply(filename);
        }
        } else if (extension == "bdf" || extension == "dat" || extension == "nas") {
            logger.info("Detected Nastran file format.");
            return load_mesh_from_nastran(filename);
        } else if (extension == "txt") {
            logger.info("Detected text file format.");
            return load_mesh_from_txt(filename);
        }
    logger.error("Unsupported mesh file format.");
    return false;
}

// Private function to parse Nastran Bulk Data Format (.bdf, .dat)
bool Mesh::load_mesh_from_nastran(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        logger.error("Failed to open file: " + filename);
        return false;
    }

    nodes.clear();
    elements.clear();

    std::string line;
    bool in_bulk_section = false;

    try {
        while (std::getline(file, line)) {
            // Ignore comment lines
            if (line.rfind('$', 0) == 0) {
                continue;
            }

            // The data we care about is between BEGIN BULK and ENDDATA
            if (line.rfind("BEGIN BULK", 0) == 0) {
                in_bulk_section = true;
                continue;
            }
            if (line.rfind("ENDDATA", 0) == 0) {
                break;
            }
            if (!in_bulk_section) {
                continue;
            }

            // Parse GRID cards for nodes (fixed small-field format)
            if (line.rfind("GRID", 0) == 0 && line.length() >= 48) {
                Node node;
                int nastran_id = std::stoi(trim(line.substr(8, 8)));
                node.id = nastran_id;

                node.coordinates.x() = std::stod(trim(line.substr(24, 8)));
                node.coordinates.y() = std::stod(trim(line.substr(32, 8)));
                node.coordinates.z() = std::stod(trim(line.substr(40, 8)));

                nodes.push_back(node);
            }
            // Parse CTETRA cards for 4-node tetrahedral elements
            else if (line.rfind("CTETRA", 0) == 0 && line.length() >= 40) {
                if (this->element_vertice_num != 4) {
                     logger.error("Attempting to load tetrahedral elements into a mesh configured for " +
                                  std::to_string(this->element_vertice_num) + " nodes per element.");
                     return false;
                }
                int element_id = std::stoi(trim(line.substr(8, 8)));
                Element element(4, element_id);
                // element.id is already set by constructor

                int g1 = std::stoi(trim(line.substr(24, 8)));
                int g2 = std::stoi(trim(line.substr(32, 8)));
                int g3 = std::stoi(trim(line.substr(40, 8)));
                int g4 = std::stoi(trim(line.substr(48, 8)));

                element.node_ids[0] = g1;
                element.node_ids[1] = g2;
                element.node_ids[2] = g3;
                element.node_ids[3] = g4;


                elements.push_back(element);
            }
        }
    } catch (const std::exception& e) {
        logger.error("An error occurred while parsing the BDF file: ", e.what());
        return false;
    }

    logger.info("Nastran BDF mesh loaded successfully.");
    logger.info("Mesh has " + std::to_string(elements.size()) + " elements and " + std::to_string(nodes.size()) + " nodes.");
    return true;
}


// --- Your existing PLY loader ---
bool Mesh::load_mesh_from_ply(const std::string &filename) {
    logger.error("PLY mesh loading not implemented.");
    return false;
    // logger.info("Parsing PLY file...");
    // std::ifstream file(filename);
    // std::string line;
    //
    // int num_nodes = 0;
    // int num_elements = 0;
    // bool header_ended = false;
    //
    // // Header parsing
    // while (std::getline(file, line)) {
    //     std::stringstream ss(line);
    //     std::string keyword;
    //     ss >> keyword;
    //
    //     if (keyword == "format") {
    //         std::string format_type;
    //         ss >> format_type;
    //         if (format_type != "ascii") {
    //             logger.error("Error: Only ascii PLY files are supported.");
    //             return false;
    //         }
    //     } else if (keyword == "element") {
    //         std::string element_type;
    //         ss >> element_type;
    //         if (element_type == "vertex") {
    //             ss >> num_nodes;
    //         } else if (element_type == "face") {
    //             ss >> num_elements;
    //         }
    //     } else if (keyword == "end_header") {
    //         header_ended = true;
    //         break;
    //     }
    // }
    //
    // if (!header_ended) {
    //     logger.error("Error: Invalid PLY file. 'end_header' not found.");
    //     return false;
    // }
    //
    // // Data parsing
    // nodes.reserve(num_nodes);
    // for (int i = 0; i < num_nodes; ++i) {
    //     if (!std::getline(file, line)) {
    //         logger.error("Error: Unexpected end of file while reading nodes.");
    //         return false;
    //     }
    //     std::stringstream ss(line);
    //     Node node;
    //     node.id = i;
    //     ss >> node.coordinates[0] >> node.coordinates[1] >> node.coordinates[2];
    //     nodes.push_back(node);
    // }
    //
    // elements.reserve(num_elements);
    // for (int i = 0; i < num_elements; ++i) {
    //     if (!std::getline(file, line)) {
    //         logger.error("Error: Unexpected end of file while reading elements.");
    //         return false;
    //     }
    //     std::stringstream ss(line);
    //     int face_verts;
    //     ss >> face_verts;
    //
    //     if (face_verts != this->element_vertice_num) {
    //         logger.error("Error: Mesh file contains elements with " + std::to_string(face_verts) +
    //                      " vertices, but solver is configured for " + std::to_string(this->element_vertice_num) + ".");
    //         return false;
    //     }
    //
    //     Element element(this->element_vertice_num, i);
    //     for (int j = 0; j < this->element_vertice_num; ++j) {
    //         ss >> element.node_ids[j];
    //     }
    //     elements.push_back(element);
    // }
    //
    // logger.info("PLY mesh loaded successfully.");
    // logger.info("Mesh has " + std::to_string(elements.size()) + " elements and " + std::to_string(nodes.size()) + " nodes.");
    // return true;
}
bool Mesh::load_mesh_from_txt(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        logger.error("Failed to open file: " + filename);
        return false;
    }

    nodes.clear();
    elements.clear();

    std::string line;
    enum class ParseState { Header, Nodes, Elements };
    ParseState state = ParseState::Header;

    int node_count = 0;
    int element_count = 0;

    try {
        while (std::getline(file, line)) {
            // Skip empty lines
            if (line.empty()) continue;

            // Simple state machine based on comment headers
            if (line.rfind("% Coordinates", 0) == 0) {
                state = ParseState::Nodes;
                continue;
            } else if (line.rfind("% Elements", 0) == 0) {
                state = ParseState::Elements;
                continue;
            }

            // Skip any other comment lines
            if (line.rfind('%', 0) == 0) {
                continue;
            }

            // Process data based on the current state
            if (state == ParseState::Nodes) {
                std::stringstream ss(line);
                Node node;
                node.id = node_count++;
                ss >> node.coordinates.x() >> node.coordinates.y() >> node.coordinates.z();
                nodes.push_back(node);
            } else if (state == ParseState::Elements) {
                if (this->element_vertice_num != 4) {
                    logger.error("Attempting to load 4-node tetrahedral elements into a mesh configured for " +
                                 std::to_string(this->element_vertice_num) + " nodes per element.");
                    return false;
                }
                std::stringstream ss(line);
                Element element(4, element_count++);
                int n1, n2, n3, n4;
                ss >> n1 >> n2 >> n3 >> n4;

                // The file format uses 1-based indexing. Convert to 0-based for our solver.
                element.node_ids[0] = n1 - 1;
                element.node_ids[1] = n2 - 1;
                element.node_ids[2] = n3 - 1;
                element.node_ids[3] = n4 - 1;

                elements.push_back(element);
            }
        }
    } catch (const std::exception& e) {
        logger.error("An error occurred while parsing the TXT file: ", e.what());
        return false;
    }

    logger.info("COMSOL TXT mesh loaded successfully.");
    logger.info("Mesh has " + std::to_string(elements.size()) + " elements and " + std::to_string(nodes.size()) + " nodes.");
    return true;
}