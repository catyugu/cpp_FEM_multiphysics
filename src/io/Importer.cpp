#include "io/Importer.hpp"
#include <core/mesh/Mesh.hpp>
#include <core/mesh/Node.hpp>
#include <core/mesh/TriElement.hpp> // For 2D meshes
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>


namespace IO {
    std::unique_ptr<Core::Mesh> Importer::read_comsol_mphtxt(const std::string &filename) {
        auto &logger = SimpleLogger::Logger::instance();
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
                        Core::Node *node1 = mesh->getNode(n1);
                        Core::Node *node2 = mesh->getNode(n2);
                        Core::Node *node3 = mesh->getNode(n3);

                        if (node1 && node2 && node3) {
                            auto *tri_elem = new Core::TriElement(element_id_counter++);
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

        logger.info("Import finished. Total: ", mesh->getNodes().size(), " nodes, ", mesh->getElements().size(),
                    " elements.");

        if (mesh->getNodes().empty() || mesh->getElements().empty()) {
            throw Exception::FileIOException(
                "The importer did not read any valid mesh data. Please check the .mphtxt file format and content.");
        }

        return mesh;
    }

    std::unique_ptr<Core::Mesh> Importer::read_gmsh_msh(const std::string &filename) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("Importing Gmsh mesh from file: ", filename);

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw Exception::FileIOException("Failed to open mesh file: " + filename);
        }

        auto mesh = std::make_unique<Core::Mesh>();
        std::string line;
        std::map<int, int> gmsh_node_id_to_internal_id;
        int internal_node_id_counter = 0;
        int element_id_counter = 0;

        while (std::getline(file, line)) {
            if (line.find("$Nodes") != std::string::npos) {
                std::getline(file, line); // Read the line with the number of nodes
                long long num_nodes;
                std::stringstream(line) >> num_nodes;

                logger.info("Reading ", num_nodes, " nodes...");
                for (long long i = 0; i < num_nodes; ++i) {
                    std::getline(file, line);
                    std::stringstream ss(line);
                    int gmsh_node_id;
                    double x, y, z;
                    ss >> gmsh_node_id >> x >> y >> z;

                    mesh->addNode(new Core::Node(internal_node_id_counter, x, y, z));
                    gmsh_node_id_to_internal_id[gmsh_node_id] = internal_node_id_counter;
                    internal_node_id_counter++;
                }
            } else if (line.find("$Elements") != std::string::npos) {
                std::getline(file, line); // Read the line with the number of elements
                long long num_elements;
                std::stringstream(line) >> num_elements;

                logger.info("Reading ", num_elements, " elements...");
                for (long long i = 0; i < num_elements; ++i) {
                    std::getline(file, line);
                    std::stringstream ss(line);
                    int elem_tag, elem_type, num_tags;
                    ss >> elem_tag >> elem_type >> num_tags;

                    // Skip tags
                    for (int j = 0; j < num_tags; ++j) {
                        int tag;
                        ss >> tag;
                    }

                    std::vector<int> node_ids;
                    int node_id;
                    while (ss >> node_id) {
                        node_ids.push_back(gmsh_node_id_to_internal_id.at(node_id));
                    }

                    Core::Element *new_elem = nullptr;
                    switch (elem_type) {
                        case 1: // 2-node line
                            if (node_ids.size() == 2) {
                                new_elem = new Core::LineElement(element_id_counter++);
                            }
                            break;
                        case 2: // 3-node triangle
                            if (node_ids.size() == 3) {
                                new_elem = new Core::TriElement(element_id_counter++);
                            }
                            break;
                        case 4: // 4-node tetrahedron
                            if (node_ids.size() == 4) {
                                new_elem = new Core::TetElement(element_id_counter++);
                            }
                            break;
                        default:
                            logger.warn("Unsupported Gmsh element type found: ", elem_type);
                            break;
                    }

                    if (new_elem) {
                        for (int id: node_ids) {
                            new_elem->addNode(mesh->getNode(id));
                        }
                        mesh->addElement(new_elem);
                    }
                }
            }
        }

        logger.info("Import finished. Total: ", mesh->getNodes().size(), " nodes, ", mesh->getElements().size(),
                    " elements.");
        if (mesh->getNodes().empty() || mesh->getElements().empty()) {
            throw Exception::FileIOException("The importer did not read any valid mesh data from the .msh file.");
        }
        return mesh;
    }

    std::vector<double> Importer::read_vtu_data(const std::string &filename, const std::string &data_array_name) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("Importing VTU data from file: ", filename, " for array: ", data_array_name);

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw Exception::FileIOException("Failed to open VTU file: " + filename);
        }

        std::vector<double> data;
        std::string line;
        bool in_correct_data_array = false;

        while (std::getline(file, line)) {
            // Find the start of the correct DataArray
            if (line.find("<DataArray") != std::string::npos && line.find("Name=\"" + data_array_name + "\"") !=
                std::string::npos) {
                in_correct_data_array = true;
                continue; // Move to the next line to start reading data
            }

            // Find the end of the DataArray
            if (in_correct_data_array && line.find("</DataArray>") != std::string::npos) {
                break; // Stop reading
            }

            // If we are in the correct block, parse the numbers
            if (in_correct_data_array) {
                std::stringstream ss(line);
                double value;
                while (ss >> value) {
                    data.push_back(value);
                }
            }
        }

        if (data.empty()) {
            logger.warn("No data was read for array '", data_array_name, "' from file '", filename, "'.");
        } else {
            logger.info("Successfully read ", data.size(), " data points for array '", data_array_name, "'.");
        }

        return data;
    }

    VtuData Importer::read_vtu_points_and_data(const std::string &filename,
                                               const std::vector<std::string> &data_array_names) {
        auto &logger = SimpleLogger::Logger::instance();
        logger.info("Importing points and data from VTU file: ", filename);

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw Exception::FileIOException("Failed to open VTU file: " + filename);
        }

        VtuData vtu_data;
        std::string line;
        bool in_points_section = false;
        std::string current_data_array_name;

        while (std::getline(file, line)) {
            // Find and read points
            if (line.find("<Points>") != std::string::npos) {
                in_points_section = true;
                std::getline(file, line); // Skip the <DataArray> line for points for now
            } else if (in_points_section && line.find("</Points>") != std::string::npos) {
                in_points_section = false;
            } else if (in_points_section && line.find("<DataArray") == std::string::npos) {
                std::stringstream ss(line);
                double x, y, z;
                ss >> x >> y >> z;
                vtu_data.points.push_back({x, y, z});
            }

            // Find and read specified data arrays
            for (const auto &name: data_array_names) {
                if (line.find("<DataArray") != std::string::npos && line.find("Name=\"" + name + "\"") !=
                    std::string::npos) {
                    current_data_array_name = name;
                    break;
                }
            }

            if (!current_data_array_name.empty() && line.find("</DataArray>") != std::string::npos) {
                current_data_array_name.clear();
            } else if (!current_data_array_name.empty() && line.find("<DataArray") == std::string::npos) {
                std::stringstream ss(line);
                double value;
                while (ss >> value) {
                    vtu_data.point_data[current_data_array_name].push_back(value);
                }
            }
        }

        logger.info("Finished reading VTU. Found ", vtu_data.points.size(), " points.");
        return vtu_data;
    }
} // namespace IO
