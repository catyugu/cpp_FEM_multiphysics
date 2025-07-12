#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/Element.hpp>
#include <core/mesh/Node.hpp>
#include "utils/Quadrature.hpp"
#include "core/sources/SourceTerm.hpp"

namespace Physics {

// ... constructor and setup are unchanged ...

void Heat1D::assemble() {
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

    std::vector<Eigen::Triplet<double>> k_triplets, m_triplets;
    auto quadrature_points = Utils::Quadrature::getLineQuadrature(element_order_);

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* line_elem = dynamic_cast<Core::LineElement*>(mesh_->getElements()[i]);
        if (line_elem) {
            double h = line_elem->getLength();

            Eigen::Matrix2d ke_local = Eigen::Matrix2d::Zero();
            Eigen::Matrix2d me_local = Eigen::Matrix2d::Zero();

            for(const auto& qp : quadrature_points) {
                double detJ = h / 2.0;
                Eigen::Matrix<double, 1, 2> B;
                B << -1/h, 1/h;
                ke_local += B.transpose() * k * B * qp.weight * detJ;

                Eigen::Matrix<double, 1, 2> N;
                N << (1.0 - qp.point(0))/2.0, (1.0 + qp.point(0))/2.0;
                me_local += N.transpose() * (rho * cp) * N * qp.weight * detJ;
            }

            auto nodes = line_elem->getNodes();
            int dofs[2] = {
                dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName()),
                dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName())
            };

            for(int r = 0; r < 2; ++r) {
                for(int c = 0; c < 2; ++c) {
                    k_triplets.emplace_back(dofs[r], dofs[c], ke_local(r,c));
                    m_triplets.emplace_back(dofs[r], dofs[c], me_local(r,c));
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());

    // THE FIX: The stabilization loop has been removed.

    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics