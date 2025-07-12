#include "physics/Magnetic2D.hpp"
#include <core/mesh/TriElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"

namespace Physics {

Magnetic2D::Magnetic2D(const Core::Material& material)
    : material_(material) {}

const char* Magnetic2D::getName() const { return "Magnetic Field 2D"; }
const char* Magnetic2D::getVariableName() const { return "MagneticPotential"; }

void Magnetic2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);

    K_.setZero();
    F_.setZero();
    U_.setZero();
}

void Magnetic2D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    const double mu = material_.getProperty("magnetic_permeability");
    const double inv_mu = 1.0 / mu;

    std::vector<Eigen::Triplet<double>> k_triplets;
    auto quadrature_points = Utils::Quadrature::getTriangleQuadrature(element_order_);

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
        if (tri_elem) {
            Eigen::Matrix3d ke_local = Eigen::Matrix3d::Zero();
            for (const auto& qp : quadrature_points) {
                auto B = tri_elem->getBMatrix();
                double detJ = tri_elem->getArea() * 2.0;
                ke_local += B.transpose() * inv_mu * B * qp.weight * detJ;
            }

            auto nodes = tri_elem->getNodes();
            int dofs[3];
            for(int i=0; i<3; ++i) {
                dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            }

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                }
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());

    const std::string my_var = getVariableName();
    for(const auto& var_name : dof_manager_->getVariableNames()) {
        if (var_name != my_var) {
            for(const auto& node : mesh_->getNodes()) {
                int dof_idx = dof_manager_->getEquationIndex(node->getId(), var_name);
                if (dof_idx != -1) K_.coeffRef(dof_idx, dof_idx) = 1.0;
            }
        }
    }
    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics