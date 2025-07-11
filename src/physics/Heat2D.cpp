#include "physics/Heat2D.hpp"
#include "core/TriElement.hpp"
#include "core/Node.hpp"
#include "utils/SimpleLogger.hpp"

namespace Physics {

// Constructor and setup are correct...
Heat2D::Heat2D(const Core::Material& material) : material_(material), k_(0.0), rho_(0.0), cp_(0.0) {}
const char* Heat2D::getName() const { return "Heat Transfer 2D"; }
const char* Heat2D::getVariableName() const { return "Temperature"; }

void Heat2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    k_ = material_.getProperty("thermal_conductivity");
    rho_ = material_.getProperty("density");
    cp_ = material_.getProperty("specific_heat");

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> k = ", k_, ", rho = ", rho_, ", cp = ", cp_);

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq,1); F_.setZero();
    U_.resize(num_eq,1); U_.setZero();
    U_prev_.resize(num_eq,1); U_prev_.setZero();
    volumetric_heat_source_.resize(mesh_->getElements().size(), 0.0);
}

void Heat2D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    M_.setZero(); // <-- FIX: Ensure Mass Matrix is also zeroed out
    F_.setZero();

    Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * k_;
    std::vector<Eigen::Triplet<double>> k_triplets;
    std::vector<Eigen::Triplet<double>> m_triplets; // <-- FIX: Add triplets for Mass Matrix

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(mesh_->getElements()[i]);
        if (tri_elem) {
            double area = tri_elem->getArea();
            auto B = tri_elem->getBMatrix();

            // Stiffness Matrix: k_e = Area * B^T * D * B
            Eigen::Matrix3d ke = area * B.transpose() * D * B;

            // --- FIX: Add Mass Matrix calculation ---
            // Mass Matrix: m_e = (rho * cp * Area / 12) * [[2,1,1],[1,2,1],[1,1,2]]
            Eigen::Matrix3d me_base;
            me_base << 2, 1, 1, 1, 2, 1, 1, 1, 2;
            Eigen::Matrix3d me = (rho_ * cp_ * area / 12.0) * me_base;

            auto nodes = tri_elem->getNodes();
            int dofs[3];
            for(int j=0; j<3; ++j) {
                dofs[j] = dof_manager_->getEquationIndex(nodes[j]->getId(), getVariableName());
            }

            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    k_triplets.emplace_back(dofs[r], dofs[c], ke(r, c));
                    m_triplets.emplace_back(dofs[r], dofs[c], me(r, c)); // <-- FIX
                }
            }

            double fe_val = volumetric_heat_source_[i] * area / 3.0;
            if (fe_val != 0.0) {
                for (int j=0; j<3; ++j) F_(dofs[j]) += fe_val;
            }
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end()); // <-- FIX: Assemble the global M matrix

    // Stabilization loop for K matrix (for coupled problems)
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

void Heat2D::setVolumetricHeatSource(const std::vector<double>& source) {
    volumetric_heat_source_ = source;
}

} // namespace Physics
