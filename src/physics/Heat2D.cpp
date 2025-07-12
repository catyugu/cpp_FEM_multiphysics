#include "physics/Heat2D.hpp"
#include <core/mesh/TriElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"
#include "core/sources/SourceTerm.hpp"


namespace Physics {

Heat2D::Heat2D(const Core::Material& material) : material_(material) {}

const char* Heat2D::getName() const { return "Heat Transfer 2D"; }
const char* Heat2D::getVariableName() const { return "Temperature"; }

void Heat2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq,1); F_.setZero();
    U_.resize(num_eq,1); U_.setZero();
    U_prev_.resize(num_eq,1); U_prev_.setZero();
}
void Heat2D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    M_.setZero();

    F_.setZero();
    for (const auto& source : source_terms_) {
        source->apply(F_, *dof_manager_, *mesh_, getVariableName());
    }

    const double k = material_.getProperty("thermal_conductivity");
    const double rho = material_.getProperty("density");
    const double cp = material_.getProperty("specific_heat");
    Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * k;

    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets;
    auto quadrature_points = Utils::Quadrature::getTriangleQuadrature(element_order_);

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(mesh_->getElements()[i]);
        if (tri_elem) {
            Eigen::Matrix3d ke_local = Eigen::Matrix3d::Zero();
            Eigen::Matrix3d me_local = Eigen::Matrix3d::Zero();
            double detJ = tri_elem->getArea() * 2.0; // Jacobian for triangle

            for (const auto& qp : quadrature_points) {
                auto B = tri_elem->getBMatrix();
                ke_local += B.transpose() * D * B * qp.weight * detJ;

                Eigen::Matrix<double, 1, 3> N;
                // Shape functions for a linear triangle in natural coordinates
                N << 1.0 - qp.point(0) - qp.point(1), qp.point(0), qp.point(1);
                me_local += N.transpose() * (rho * cp) * N * qp.weight * detJ;
            }

            auto nodes = tri_elem->getNodes();
            int dofs[3];
            for(int j=0; j<3; ++j) {
                dofs[j] = dof_manager_->getEquationIndex(nodes[j]->getId(), getVariableName());
            }

            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    k_triplets.emplace_back(dofs[r], dofs[c], ke_local(r, c));
                    m_triplets.emplace_back(dofs[r], dofs[c], me_local(r, c));
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());

    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics