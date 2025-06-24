#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Element.hpp"
#include "core/Node.hpp"

namespace Physics {

Heat1D::Heat1D(const Core::Material& material)
    : material_(material), k_(0.0), rho_(0.0), cp_(0.0) {}

const char* Heat1D::getName() const { return "Heat Transfer 1D"; }
const char* Heat1D::getVariableName() const { return "Temperature"; }


void Heat1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
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
    F_.resize(num_eq); F_.setZero();
    U_.resize(num_eq); U_.setZero();
    U_prev_.resize(num_eq); U_prev_.setZero();
    volumetric_heat_source_.resize(mesh_->getElements().size(), 0.0);
}

void Heat1D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    M_.setZero();
    F_.setZero();

    std::vector<Eigen::Triplet<double>> k_triplets, m_triplets;

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* line_elem = dynamic_cast<Core::LineElement*>(mesh_->getElements()[i]);
        if (line_elem) {
            double h = line_elem->getLength();
            double ke_val = k_ / h;

            double me_val_diag = (rho_ * cp_ * h / 6.0) * 2.0;
            double me_val_offdiag = (rho_ * cp_ * h / 6.0) * 1.0;

            auto nodes = line_elem->getNodes();
            int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
            int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());

            k_triplets.emplace_back(dof_i, dof_i, ke_val);
            k_triplets.emplace_back(dof_i, dof_j, -ke_val);
            k_triplets.emplace_back(dof_j, dof_i, -ke_val);
            k_triplets.emplace_back(dof_j, dof_j, ke_val);

            m_triplets.emplace_back(dof_i, dof_i, me_val_diag);
            m_triplets.emplace_back(dof_i, dof_j, me_val_offdiag);
            m_triplets.emplace_back(dof_j, dof_i, me_val_offdiag);
            m_triplets.emplace_back(dof_j, dof_j, me_val_diag);

            double fe_val = volumetric_heat_source_[i] * h / 2.0;
            F_(dof_i) += fe_val;
            F_(dof_j) += fe_val;
        }
    }
    K_.setFromTriplets(k_triplets.begin(), k_triplets.end());
    M_.setFromTriplets(m_triplets.begin(), m_triplets.end());

    // --- FIX: Add stabilization loop for other physics' DOFs ---
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

void Heat1D::setVolumetricHeatSource(const std::vector<double>& source) {
    if (source.size() != mesh_->getElements().size()) {
        SimpleLogger::Logger::instance().error("Heat source vector size mismatch in Heat1D.");
        return;
    }
    volumetric_heat_source_ = source;
}

} // namespace Physics
