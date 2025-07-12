#include "physics/Heat3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"

namespace Physics {

Heat3D::Heat3D(const Core::Material& material) : material_(material), k_(0.0) {}

const char* Heat3D::getName() const { return "Heat Transfer 3D"; }
const char* Heat3D::getVariableName() const { return "Temperature"; }

void Heat3D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    k_ = material_.getProperty("thermal_conductivity");

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> Thermal Conductivity (k): ", k_);

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);
    F_.setZero();
}

void Heat3D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    // Material matrix for 3D isotropic heat conduction
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * k_;

    std::vector<Eigen::Triplet<double>> triplet_list;

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            double volume = tet_elem->getVolume();
            auto B = tet_elem->getBMatrix();

            // Element stiffness matrix: k_e = Volume * B^T * D * B
            Eigen::Matrix4d ke = volume * B.transpose() * D * B;

            auto nodes = tet_elem->getNodes();
            int dofs[4];
            for(int i=0; i<4; ++i) {
                dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            }

            // Add element matrix to global matrix triplets
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    triplet_list.emplace_back(dofs[i], dofs[j], ke(i, j));
                }
            }
        }
    }
    K_.setFromTriplets(triplet_list.begin(), triplet_list.end());
    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics
