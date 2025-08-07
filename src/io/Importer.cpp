#include "io/Importer.hpp"
#include <core/mesh/Mesh.hpp>
#include <core/mesh/Node.hpp>
#include <core/mesh/TriElement.hpp>
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Exceptions.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <core/mesh/LineElement.hpp>


namespace IO {
    std::unique_ptr<Core::Mesh> Importer::read_comsol_mphtxt(const std::string &filename) {
        auto &logger = Utils::Logger::instance();
        logger.info("Importing COMSOL mesh from file: ", filename);

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw Exception::FileIOException("Failed to open mesh file: " + filename);
        }

        auto mesh = std::make_unique<Core::Mesh>();
        int element_id_counter = 0;
        int sdim = 0; // 存储空间维度
        std::string line;

        file.clear();
        file.seekg(0, std::ios::beg);

        // --- 预扫描: 找到空间维度 sdim ---
        while (std::getline(file, line)) {
            if (line.find("# sdim") != std::string::npos) {
                std::stringstream ss(line);
                ss >> sdim;
                logger.info("Detected spatial dimension (sdim) = ", sdim);
                break;
            }
        }
        if (sdim == 0) {
            throw Exception::FileIOException("Could not determine spatial dimension (sdim) from mesh file.");
        }


        // --- Pass 1: 读取所有顶点/节点 ---
        file.clear();
        file.seekg(0, std::ios::beg);
        logger.info("Pass 1: Reading mesh vertices...");
        while (std::getline(file, line)) {
            if (line.find("number of mesh vertices") != std::string::npos) {
                int num_vertices = 0;
                std::stringstream ss(line);
                ss >> num_vertices;

                while (std::getline(file, line) && line.find("# Mesh vertex coordinates") == std::string::npos) {}

                for (int i = 0; i < num_vertices && std::getline(file, line); ++i) {
                    std::stringstream data_ss(line);
                    double x = 0.0, y = 0.0, z = 0.0;
                    data_ss >> x >> y >> z;
                    int node_id = static_cast<int>(mesh->getNodes().size());
                    mesh->addNode(new Core::Node(node_id, x, y, z));
                }
                logger.info("Finished reading ", mesh->getNodes().size(), " vertices.");
                break;
            }
        }

        // --- Pass 2: 只读取与问题维度相符的单元 ---
        file.clear();
        file.seekg(0, std::ios::beg);
        logger.info("Pass 2: Reading mesh elements for ", sdim, "D problem...");
        int loaded_elements = 0;
        while (std::getline(file, line)) {
            // 如果是2D问题，只读取三角形单元
            if (sdim == 2 && line.find("tri # type name") != std::string::npos) {
                while (std::getline(file, line) && line.find("# number of elements") == std::string::npos) {}
                int num_elements = 0;
                std::stringstream ss(line);
                ss >> num_elements;
                while (std::getline(file, line) && line.find("# Elements") == std::string::npos) {}
                for (int i = 0; i < num_elements && std::getline(file, line); ++i) {
                    std::stringstream data_ss(line);
                    int n1, n2, n3;
                    if (data_ss >> n1 >> n2 >> n3) {
                        auto *tri_elem = new Core::TriElement(element_id_counter++);
                        tri_elem->addNode(mesh->getNode(n1));
                        tri_elem->addNode(mesh->getNode(n2));
                        tri_elem->addNode(mesh->getNode(n3));
                        tri_elem->update_geometry();
                        mesh->addElement(tri_elem);
                        loaded_elements++;
                    }
                }
            }
            // 如果是3D问题，只读取四面体单元
            else if (sdim == 3 && line.find("tet # type name") != std::string::npos) {
                while (std::getline(file, line) && line.find("# number of elements") == std::string::npos) {}
                int num_elements = 0;
                std::stringstream ss(line);
                ss >> num_elements;
                while (std::getline(file, line) && line.find("# Elements") == std::string::npos) {}
                for (int i = 0; i < num_elements && std::getline(file, line); ++i) {
                    std::stringstream data_ss(line);
                    int n1, n2, n3, n4;
                    if (data_ss >> n1 >> n2 >> n3 >> n4) {
                        auto *tet_elem = new Core::TetElement(element_id_counter++);
                         std::vector<Core::Node*> current_nodes_order = {
                            mesh->getNode(n1), mesh->getNode(n2), mesh->getNode(n3), mesh->getNode(n4)
                        };
                        for (Core::Node* node_ptr : current_nodes_order) {
                            tet_elem->addNode(node_ptr);
                        }
                        tet_elem->update_geometry();
                        double volume = tet_elem->getVolume();

                        if (std::abs(volume) < 1e-12) {
                            logger.warn("Degenerate (zero volume) element ", tet_elem->getId(), " skipped.");
                            delete tet_elem;
                            continue;
                        }

                        if (volume < 0) {
                             std::swap(current_nodes_order[2], current_nodes_order[3]);
                             tet_elem->set_nodes_internal(current_nodes_order);
                             tet_elem->update_geometry();
                             if (tet_elem->getVolume() < 0) {
                                logger.error("Inverted element ", tet_elem->getId(), " could not be fixed. Skipping.");
                                delete tet_elem;
                                continue;
                             }
                        }
                        mesh->addElement(tet_elem);
                        loaded_elements++;
                    }
                }
            }
        }

        logger.info("Import finished. Loaded ", loaded_elements, " dimensionally-correct elements.");
        logger.info("Total: ", mesh->getNodes().size(), " nodes, ", mesh->getElements().size(), " elements.");

        if (mesh->getNodes().empty() || mesh->getElements().empty()) {
            throw Exception::FileIOException(
                "The importer did not read any valid mesh data. Please check the .mphtxt file format and content.");
        }

        return mesh;
    }

    // ... (rest of the file remains unchanged) ...
    std::unique_ptr<Core::Mesh> Importer::read_gmsh_msh(const std::string &filename) {
        auto &logger = Utils::Logger::instance();
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
                        new_elem->update_geometry();
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
        auto &logger = Utils::Logger::instance();
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
        auto &logger = Utils::Logger::instance();
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
}