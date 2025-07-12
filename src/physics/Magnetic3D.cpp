#include "physics/Magnetic3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"

namespace Physics {

Magnetic3D::Magnetic3D(const Core::Material& material) : material_(material) {}

const char* Magnetic3D::getName() const { return "Magnetic Field 3D"; }
const char* Magnetic3D::getVariableName() const { return "MagneticPotential"; }

void Magnetic3D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);
    F_.setZero();
}

void Magnetic3D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    const double mu = material_.getProperty("magnetic_permeability");
    const double inv_mu = 1.0 / mu;
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * inv_mu;

    std::vector<Eigen::Triplet<double>> triplet_list;
    auto quadrature_points = Utils::Quadrature::getTetrahedronQuadrature(element_order_);

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            Eigen::Matrix4d ke_local = Eigen::Matrix4d::Zero();
            for (const auto& qp : quadrature_points) {
                auto B = tet_elem->getBMatrix();
                double detJ = tet_elem->getVolume() * 6.0;
                ke_local += B.transpose() * D * B * qp.weight * detJ;
            }

            auto nodes = tet_elem->getNodes();
            int dofs[4];
            for(int i=0; i<4; ++i) {
                dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    triplet_list.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                }
            }
        }
    }
    K_.setFromTriplets(triplet_list.begin(), triplet_list.end());
    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics