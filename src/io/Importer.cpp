#include "io/Importer.hpp"
#include <core/mesh/Mesh.hpp>
#include <core/mesh/Node.hpp>
#include <core/mesh/TriElement.hpp> // For 2D meshes
#include <core/mesh/TetElement.hpp> // For 3D meshes
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
        int element_id_counter = 0;
        std::string line;

        // --- Pass 1: Read all Vertices/Nodes ---
        file.clear();
        file.seekg(0, std::ios::beg);
        logger.info("Pass 1: Reading mesh vertices...");
        while (std::getline(file, line)) {
            if (line.find("number of mesh vertices") != std::string::npos) {
                int num_vertices = 0;
                std::stringstream ss(line);
                ss >> num_vertices;

                while (std::getline(file, line) && line.find("# Mesh vertex coordinates") == std::string::npos);

                for (int i = 0; i < num_vertices && std::getline(file, line); ++i) {
                    std::stringstream data_ss(line);
                    double x, y, z = 0.0;
                    data_ss >> x >> y >> z;
                    int node_id = mesh->getNodes().size();
                    mesh->addNode(new Core::Node(node_id, x, y, z));
                }
                logger.info("Finished reading ", mesh->getNodes().size(), " vertices.");
                break;
            }
        }

        // --- Pass 2: Read all Triangle Elements (if any) ---
        file.clear();
        file.seekg(0, std::ios::beg);
        logger.info("Pass 2: Reading triangle elements...");
        while (std::getline(file, line)) {
            if (line.find("tri # type name") != std::string::npos) {
                while (std::getline(file, line) && line.find("# number of elements") == std::string::npos);
                int num_elements = 0;
                std::stringstream ss(line);
                ss >> num_elements;

                while (std::getline(file, line) && line.find("# Elements") == std::string::npos);

                for (int i = 0; i < num_elements && std::getline(file, line); ++i) {
                    std::stringstream data_ss(line);
                    int n1, n2, n3;
                    if (data_ss >> n1 >> n2 >> n3) {
                        auto *tri_elem = new Core::TriElement(element_id_counter++);
                        tri_elem->addNode(mesh->getNode(n1));
                        tri_elem->addNode(mesh->getNode(n2));
                        tri_elem->addNode(mesh->getNode(n3));
                        mesh->addElement(tri_elem);
                    }
                }
                break;
            }
        }

        // --- Pass 3: Read all Tetrahedron Elements (if any) ---
        file.clear();
        file.seekg(0, std::ios::beg);
        logger.info("Pass 3: Reading tetrahedron elements...");
        while (std::getline(file, line)) {
            if (line.find("tet # type name") != std::string::npos) {
                while (std::getline(file, line) && line.find("# number of elements") == std::string::npos);
                int num_elements = 0;
                std::stringstream ss(line);
                ss >> num_elements;

                while (std::getline(file, line) && line.find("# Elements") == std::string::npos);

                for (int i = 0; i < num_elements && std::getline(file, line); ++i) {
                    std::stringstream data_ss(line);
                    int n1, n2, n3, n4;
                    if (data_ss >> n1 >> n2 >> n3 >> n4) {
                        auto *tet_elem = new Core::TetElement(element_id_counter++);
                        tet_elem->addNode(mesh->getNode(n1));
                        tet_elem->addNode(mesh->getNode(n2));
                        tet_elem->addNode(mesh->getNode(n3));
                        tet_elem->addNode(mesh->getNode(n4));
                        mesh->addElement(tet_elem);
                    }
                }
                break;
            }
        }

        logger.info("Import finished. Total: ", mesh->getNodes().size(), " nodes, ", mesh->getElements().size(), " elements.");

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
                std::getline(file, line);
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
                std::getline(file, line);
                long long num_elements;
                std::stringstream(line) >> num_elements;

                logger.info("Reading ", num_elements, " elements...");
                for (long long i = 0; i < num_elements; ++i) {
                    std::getline(file, line);
                    std::stringstream ss(line);
                    int elem_tag, elem_type, num_tags;
                    ss >> elem_tag >> elem_type >> num_tags;

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
                        case 1:
                            if (node_ids.size() == 2) {
                                new_elem = new Core::LineElement(element_id_counter++);
                            }
                            break;
                        case 2:
                            if (node_ids.size() == 3) {
                                new_elem = new Core::TriElement(element_id_counter++);
                            }
                            break;
                        case 4:
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
            if (line.find("<DataArray") != std::string::npos && line.find("Name=\"" + data_array_name + "\"") !=
                std::string::npos) {
                in_correct_data_array = true;
                continue;
            }

            if (in_correct_data_array && line.find("</DataArray>") != std::string::npos) {
                break;
            }

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
            if (line.find("<Points>") != std::string::npos) {
                in_points_section = true;
                std::getline(file, line);
            } else if (in_points_section && line.find("</Points>") != std::string::npos) {
                in_points_section = false;
            } else if (in_points_section && line.find("<DataArray") == std::string::npos) {
                std::stringstream ss(line);
                double x, y, z;
                ss >> x >> y >> z;
                vtu_data.points.push_back({x, y, z});
            }

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
