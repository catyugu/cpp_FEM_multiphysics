#include <HeatField.h>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// Constructor initializes the heat field solver with a given mesh
HeatField::HeatField(std::unique_ptr<Mesh> mesh_ptr) : mesh(std::move(mesh_ptr)) {
    if (!mesh) {
        throw std::invalid_argument("Mesh pointer cannot be null.");
    }
    int num_nodes = static_cast<int>(mesh->nodes.size());
    if (num_nodes == 0) {
        throw std::runtime_error("Mesh contains no nodes.");
    }
    logger.info("Initializing heat field solver with mesh containing %d nodes.", num_nodes);
    init_condition.resize(num_nodes);
    border_condition.resize(num_nodes);
    stiffness_matrix.resize(num_nodes, num_nodes);
    mass_matrix.resize(num_nodes, num_nodes);
    F.resize(num_nodes);
    heat.resize(num_nodes);
    is_border_node.resize(num_nodes, false);

    stiffness_matrix.setZero();
    mass_matrix.setZero();
    F.setZero();
    are_matrices_assembled = false;

    // identify_border_nodes();
    logger.info("HeatField initialized for a mesh with " + std::to_string(num_nodes) + " nodes.");
}

// A simple heuristic to identify border nodes: find nodes on the mesh's bounding box.
void HeatField::identify_border_nodes() {
    if (mesh->nodes.empty()) return;

    Eigen::Vector3d min_coords = mesh->nodes[0].coordinates;
    Eigen::Vector3d max_coords = mesh->nodes[0].coordinates;

    for (const auto& node : mesh->nodes) {
        min_coords = min_coords.cwiseMin(node.coordinates);
        max_coords = max_coords.cwiseMax(node.coordinates);
    }
    logger.info("Mesh Bounding Box: (" + std::to_string(min_coords.x()) + ", " + std::to_string(min_coords.y()) + ") to (" +
                std::to_string(max_coords.x()) + ", " + std::to_string(max_coords.y()) + ")");


    constexpr double tolerance = 1e-6;
    for (size_t i = 0; i < mesh->nodes.size(); ++i) {
        const auto& coords = mesh->nodes[i].coordinates;
        if (std::abs(coords.x() - min_coords.x()) < tolerance ||
            std::abs(coords.x() - max_coords.x()) < tolerance ||
            std::abs(coords.y() - min_coords.y()) < tolerance ||
            std::abs(coords.y() - max_coords.y()) < tolerance) {
            is_border_node[i] = true;
        }
    }
}

// Assembles the global stiffness and mass matrices, and the force vector.
void HeatField::assemble_matrices() {
    logger.info("Assembling 3D stiffness and mass matrices...");
    stiffness_matrix.setZero();
    mass_matrix.setZero();
    F.setZero();

    // This implementation is for 4-node linear tetrahedral elements.
    if (mesh->element_vertice_num != 4) {
        throw std::runtime_error("3D FEM assembly is only implemented for 4-node tetrahedral elements.");
    }

    for (const auto& element : mesh->elements) {
        // 1. Get nodal coordinates for the tetrahedron
        Eigen::Vector3d p1 = mesh->nodes[element.node_ids[0]].coordinates;
        Eigen::Vector3d p2 = mesh->nodes[element.node_ids[1]].coordinates;
        Eigen::Vector3d p3 = mesh->nodes[element.node_ids[2]].coordinates;
        Eigen::Vector3d p4 = mesh->nodes[element.node_ids[3]].coordinates;

        // 2. Calculate the volume of the tetrahedron
        double volume = std::abs((p1 - p4).dot((p2 - p4).cross(p3 - p4))) / 6.0;
        if (volume < 1e-12) {
             // Skip degenerate elements with zero or negative volume
            continue;
        }

        // 3. Calculate shape function gradients (B matrix)
        // For a linear tetrahedron, the gradients are constant within the element.
        Eigen::Matrix4d M_coords;
        M_coords << 1, p1.x(), p1.y(), p1.z(),
                    1, p2.x(), p2.y(), p2.z(),
                    1, p3.x(), p3.y(), p3.z(),
                    1, p4.x(), p4.y(), p4.z();

        Eigen::Matrix4d M_inv = M_coords.inverse();
        Eigen::Matrix<double, 3, 4> B = M_inv.bottomRows<3>();

        // 4. Calculate local stiffness matrix: k_e = volume * k * B^T * B
        Eigen::Matrix4d k_e = params.k * volume * (B.transpose() * B);

        // 5. Calculate local mass matrix (consistent formulation)
        Eigen::Matrix4d m_e;
        m_e.setConstant(1.0); // Set all off-diagonal elements to 1
        m_e.diagonal().setConstant(2.0); // Set diagonal elements to 2
        m_e *= (params.rou * params.Cp * volume / 20.0);

        // 6. Calculate local force vector if source term exists
        Eigen::Vector4d f_e = Eigen::Vector4d::Zero();
        if (force_heat) {
            Eigen::Vector3d centroid = (p1 + p2 + p3 + p4) / 4.0;
            double source_at_centroid = force_heat(centroid);
            f_e << 1, 1, 1, 1;
            f_e *= source_at_centroid * volume / 4.0;
        }

        // 7. Assemble local matrices into global matrices
        for (int i = 0; i < 4; ++i) {
            F(element.node_ids[i]) += f_e(i);
            for (int j = 0; j < 4; ++j) {
                stiffness_matrix.coeffRef(element.node_ids[i], element.node_ids[j]) += k_e(i, j);
                mass_matrix.coeffRef(element.node_ids[i], element.node_ids[j]) += m_e(i, j);
            }
        }
    }
    are_matrices_assembled = true;
    logger.info("3D Matrix assembly complete.");
}

void HeatField::apply_init_condition(const std::function<double(const Eigen::Vector3d& coord)>& init_func) {
    for (size_t i = 0; i < mesh->nodes.size(); ++i) {
        heat(i) = init_func(mesh->nodes[i].coordinates);
    }
    logger.info("Initial condition applied.");
}

void HeatField::apply_border_condition(const std::function<double(const Eigen::Vector3d& coord)>& border_func) {
    for (size_t i = 0; i < mesh->nodes.size(); ++i) {
        if (is_border_node[i]) {
            border_condition(i) = border_func(mesh->nodes[i].coordinates);
        }
    }
    logger.info("Border condition values stored.");
}

void HeatField::apply_extern_field(const std::function<double(const Eigen::Vector3d& coord)>& source_func) {
    this->force_heat = source_func;
    are_matrices_assembled = false;
    logger.info("External heat source field applied.");
}

void HeatField::step_forward(double time_step) {
    if (time_step <= 0) {
        throw std::invalid_argument("Time step must be positive.");
    }
    if (!are_matrices_assembled) {
        assemble_matrices();
    }

    // Using Backward Euler: (M/dt + K) * T_next = M/dt * T_prev + F
    Eigen::MatrixXd A = mass_matrix / time_step + stiffness_matrix;
    Eigen::VectorXd b = (mass_matrix * heat) / time_step + F;

    // ====================================================================
    // CORRECTED BOUNDARY CONDITION APPLICATION
    // ====================================================================

    // Step 1: Adjust the right-hand side vector 'b' for all non-boundary nodes.
    // This "lifts" the influence of the fixed boundary temperatures over to the RHS.
    for (int i = 0; i < mesh->nodes.size(); ++i) {
        if (is_border_node[i]) {
            double boundary_value = border_condition(i);
            for (int j = 0; j < mesh->nodes.size(); ++j) {
                // For every equation 'j' that is NOT a boundary equation,
                // subtract the influence of the known boundary node 'i'.
                if (!is_border_node[j]) {
                    b(j) -= A.coeff(j, i) * boundary_value;
                }
            }
        }
    }

    // Step 2: Now that the system is adjusted, modify the matrix and RHS vector
    // to enforce the boundary condition directly.
    for (int i = 0; i < mesh->nodes.size(); ++i) {
        if (is_border_node[i]) {
            // Zero out the row and column corresponding to the boundary node
            A.row(i).setZero();
            A.col(i).setZero();

            // Place a 1 on the diagonal and the boundary value in the RHS vector
            A.coeffRef(i, i) = 1.0;
            b(i) = border_condition(i);
        }
    }

    // ====================================================================

    // Solve the linear system A * T_next = b
    // PartialPivLU is still the correct solver because the final matrix is not symmetric.
    Eigen::SparseMatrix<double> A_sparse = A.sparseView();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A_sparse);
    if(solver.info()!=Eigen::Success) {
        logger.error("LU decomposition failed.");
        return;
    }
    heat = solver.solve(b);
    if(solver.info()!=Eigen::Success) {
        logger.error("Linear system solving failed.");
        return;
    }

    logger.info("Time step " + std::to_string(time_step) + " completed.");
}

void HeatField::write_vtk(const std::string& filename) const {
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
    for (int i = 0; i < heat.size(); ++i) {
        vtk_file << heat(i) << "\n";
    }

    logger.info("Successfully wrote VTK file.");
}