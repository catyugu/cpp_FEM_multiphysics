#include "Analysis.h"
#include <fstream>
#include <iomanip>
#include <Eigen/SparseLU>

Analysis::Analysis(std::unique_ptr<Mesh> mesh_ptr) : mesh(std::move(mesh_ptr)) {
    if (!mesh) {
        throw std::invalid_argument("Mesh pointer cannot be null.");
    }
    int num_nodes = static_cast<int>(mesh->nodes.size());
    
    solution.resize(num_nodes);
    boundary_values.resize(num_nodes);
    is_border_node.resize(num_nodes, false);
    global_stiffness.resize(num_nodes, num_nodes);
    global_mass.resize(num_nodes, num_nodes);
    //
    // identify_border_nodes();
    logger.info("Analysis class initialized for a mesh with " + std::to_string(num_nodes) + " nodes.");
}

// Creates the specific physics element objects for the mesh.
void Analysis::initialize_physics() {
    logger.info("Initializing physics elements...");
    elements.clear();
    elements.reserve(mesh->elements.size());
    for (size_t i = 0; i < mesh->elements.size(); ++i) {
        // For now, we are only creating HeatTetrahedron elements.
        // This could be extended to create different elements for different mesh regions.
        elements.push_back(std::make_unique<HeatTetrahedron>());
    }
}

// Generic assembly loop, decoupled from the physics.
void Analysis::assemble_global_matrices() {
    logger.info("Assembling global matrices (Sparse)...");
    
    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets;
    
    Eigen::MatrixXd k_local(mesh->element_vertice_num, mesh->element_vertice_num);
    Eigen::MatrixXd m_local(mesh->element_vertice_num, mesh->element_vertice_num);
    
    for (size_t i = 0; i < mesh->elements.size(); ++i) {
        // Get pointers to the nodes of the current element
        std::vector<const Node*> element_nodes(mesh->element_vertice_num);
        for(int j=0; j<mesh->element_vertice_num; ++j) {
            element_nodes[j] = &mesh->nodes[mesh->elements[i].node_ids[j]];
        }
        
        // Polymorphically call the correct method to get local matrices
        elements[i]->calculate_local_matrices(element_nodes.data(), this->params, k_local, m_local);

        // Add local contributions to triplet lists
        for (int row = 0; row < mesh->element_vertice_num; ++row) {
            for (int col = 0; col < mesh->element_vertice_num; ++col) {
                int global_row = mesh->elements[i].node_ids[row];
                int global_col = mesh->elements[i].node_ids[col];
                k_triplets.emplace_back(global_row, global_col, k_local(row, col));
                m_triplets.emplace_back(global_row, global_col, m_local(row, col));
            }
        }
    }

    global_stiffness.setFromTriplets(k_triplets.begin(), k_triplets.end());
    global_mass.setFromTriplets(m_triplets.begin(), m_triplets.end());
    logger.info("Global matrix assembly complete.");
}

void Analysis::run_simulation(int num_steps, double time_step) {
    logger.info("Running simulation for " + std::to_string(num_steps) + " steps with dt=" + std::to_string(time_step));
    assemble_global_matrices();

    for (int i = 0; i < num_steps; ++i) {
        solve_step(time_step);
        logger.info("Completed step " + std::to_string(i + 1) + "/" + std::to_string(num_steps));
    }
}

void Analysis::solve_step(double time_step) {
    // This is the robust solver logic from our previous discussion
    Eigen::SparseMatrix<double> A_initial = global_mass / time_step + global_stiffness;
    Eigen::VectorXd b = (global_mass * solution) / time_step;

    std::vector<Eigen::Triplet<double>> final_triplets;
    final_triplets.reserve(A_initial.nonZeros());

    for (int k = 0; k < A_initial.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_initial, k); it; ++it) {
            if (!is_border_node[it.row()]) {
                if (is_border_node[it.col()]) {
                    b(it.row()) -= it.value() * boundary_values(it.col());
                } else {
                    final_triplets.emplace_back(it.row(), it.col(), it.value());
                }
            }
        }
    }

    for (size_t i = 0; i < is_border_node.size(); ++i) {
        if (is_border_node[i]) {
            final_triplets.emplace_back(i, i, 1.0);
            b(i) = boundary_values(i);
        }
    }

    Eigen::SparseMatrix<double> A_final(A_initial.rows(), A_initial.cols());
    A_final.setFromTriplets(final_triplets.begin(), final_triplets.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A_final);
    solver.factorize(A_final);
    if (solver.info() != Eigen::Success) {
        logger.error("Sparse LU decomposition failed.");
        return;
    }
    solution = solver.solve(b);
}

// Other helper functions (these are moved from the old HeatField class)
void Analysis::assemble_initial_condition(const std::function<double(const Eigen::Vector3d&)>& init_func) {
    for (size_t i = 0; i < mesh->nodes.size(); ++i) {
        solution(i) = init_func(mesh->nodes[i].coordinates);
    }
}
void Analysis::assemble_boundary_condition(const std::function<double(const Eigen::Vector3d&)>& border_func) {
    for (size_t i = 0; i < mesh->nodes.size(); ++i) {
        if (is_border_node[i]) {
            boundary_values(i) = border_func(mesh->nodes[i].coordinates);
        }
    }
}
void Analysis::write_vtk(const std::string& filename) const {
    logger.info("Writing results to " + filename);
    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        logger.error("Failed to open VTK file for writing: " + filename);
        return;
    }

    vtk_file << std::fixed << std::setprecision(8);

    // Header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Heat Field Simulation Result\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    // Points
    vtk_file << "POINTS " << mesh->nodes.size() << " double\n";
    for (const auto& node : mesh->nodes) {
        vtk_file << node.coordinates.x() << " " << node.coordinates.y() << " " << node.coordinates.z() << "\n";
    }

    // Cells
    const int vertices_per_element = mesh->element_vertice_num;
    vtk_file << "\nCELLS " << mesh->elements.size() << " " << mesh->elements.size() * (vertices_per_element + 1) << "\n";
    for (const auto& element : mesh->elements) {
        vtk_file << vertices_per_element;
        for (int i = 0; i < vertices_per_element; ++i) {
            vtk_file << " " << element.node_ids[i];
        }
        vtk_file << "\n";
    }

    // Cell Types (VTK_TRIANGLE = 5)
    vtk_file << "\nCELL_TYPES " << mesh->elements.size() << "\n";
    for (size_t i = 0; i < mesh->elements.size(); ++i) {
        vtk_file << "5\n";
    }

    // Point Data (Temperature)
    vtk_file << "\nPOINT_DATA " << mesh->nodes.size() << "\n";
    vtk_file << "SCALARS temperature double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < solution.size(); ++i) {
        vtk_file << solution(i) << "\n";
    }

    logger.info("Successfully wrote VTK file.");
}
// In Analysis.cpp

// Identifies border nodes by checking if they lie on the mesh's bounding box.
void Analysis::identify_border_nodes() {
    if (mesh->nodes.empty()) return;

    // Find the min/max coordinates to define a bounding box
    Eigen::Vector3d min_coords = mesh->nodes[0].coordinates;
    Eigen::Vector3d max_coords = mesh->nodes[0].coordinates;

    for (const auto& node : mesh->nodes) {
        min_coords = min_coords.cwiseMin(node.coordinates);
        max_coords = max_coords.cwiseMax(node.coordinates);
    }

    logger.info("Mesh Bounding Box Min: (" + std::to_string(min_coords.x()) + ", " + std::to_string(min_coords.y()) + ", " + std::to_string(min_coords.z()) + ")");
    logger.info("Mesh Bounding Box Max: (" + std::to_string(max_coords.x()) + ", " + std::to_string(max_coords.y()) + ", " + std::to_string(max_coords.z()) + ")");

    const double tolerance = 1e-6;
    int border_node_count = 0;
    for (size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& coords = mesh->nodes[i].coordinates;
        // Check if the node is on any of the 6 faces of the bounding box
        if (std::abs(coords.x() - min_coords.x()) < tolerance ||
            std::abs(coords.x() - max_coords.x()) < tolerance ||
            std::abs(coords.y() - min_coords.y()) < tolerance ||
            std::abs(coords.y() - max_coords.y()) < tolerance ||
            std::abs(coords.z() - min_coords.z()) < tolerance ||
            std::abs(coords.z() - max_coords.z()) < tolerance)
        {
            is_border_node[i] = true;
            border_node_count++;
        }
    }
    logger.info("Identified " + std::to_string(border_node_count) + " border nodes.");
}