#include "physics/Heat3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"
#include "utils/Exceptions.hpp"

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
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);
    F_.setZero();
}

void Heat3D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    M_.setZero();
    // F_ is handled by applySources()

    Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * k_;
    const double rho_cp = material_.getProperty("density") * material_.getProperty("specific_heat");

    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets;
    auto quadrature_points = Utils::Quadrature::getTetrahedronQuadrature(element_order_);

    int valid_elements_found = 0;
    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            valid_elements_found++;
            Eigen::Matrix4d ke_local = Eigen::Matrix4d::Zero();
            Eigen::Matrix4d me_local = Eigen::Matrix4d::Zero();
            double detJ = tet_elem->getVolume() * 6.0;

            for(const auto& qp : quadrature_points) {
                auto B = tet_elem->getBMatrix();
                ke_local += B.transpose() * D * B * qp.weight * detJ;

                Eigen::Matrix<double, 1, 4> N;
                N << 1.0 - qp.point(0) - qp.point(1) - qp.point(2), qp.point(0), qp.point(1), qp.point(2);
                me_local += N.transpose() * rho_cp * N * qp.weight * detJ;
            }

            auto nodes = tet_elem->getNodes();
            int dofs[4];
            for(int i=0; i<4; ++i) {
                dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                        m_triplets.emplace_back(dofs[i], dofs[j], me_local(i, j));
                    }
                }
            }
        }
    }

    if (valid_elements_found == 0) {
        throw Exception::ConfigurationException(
            "Assembly failed for " + std::string(getName()) +
            ": No valid elements (TetElement) were found in the mesh. Check if the mesh is a 3D volume mesh."
        );
    }

    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());

    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics
