#include "physics/EMag1D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Element.hpp"
#include "core/Node.hpp"
#include <vector>

namespace Physics {

EMag1D::EMag1D(const Core::Material& material) : material_(material), sigma_(0.0) {}

const char* EMag1D::getName() const { return "Electromagnetics 1D"; }
const char* EMag1D::getVariableName() const { return "Voltage"; }


void EMag1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;

    sigma_ = material_.getProperty("electrical_conductivity");

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> Electrical Conductivity (sigma): ", sigma_);

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq);
    U_.resize(num_eq);
}

void EMag1D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    std::vector<Eigen::Triplet<double>> triplet_list;

    for (const auto& elem_ptr : mesh_->getElements()) {
        Core::LineElement* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
        if (line_elem) {
            double h = line_elem->getLength();
            double ke_val = sigma_ / h;

            auto nodes = line_elem->getNodes();
            int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
            int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());

            triplet_list.emplace_back(dof_i, dof_i, ke_val);
            triplet_list.emplace_back(dof_i, dof_j, -ke_val);
            triplet_list.emplace_back(dof_j, dof_i, -ke_val);
            triplet_list.emplace_back(dof_j, dof_j, ke_val);
        }
    }
    K_.setFromTriplets(triplet_list.begin(), triplet_list.end());

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

std::vector<double> EMag1D::calculateJouleHeat() const {
    std::vector<double> joule_heat(mesh_->getElements().size(), 0.0);
    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        const auto& elem_ptr = mesh_->getElements()[i];
        Core::LineElement* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
        if (line_elem) {
            double h = line_elem->getLength();
            auto nodes = line_elem->getNodes();
            int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
            int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());

            double V_i = U_(dof_i);
            double V_j = U_(dof_j);
            double E = std::abs(V_i - V_j) / h;
            joule_heat[i] = sigma_ * E * E;
        }
    }
    return joule_heat;
}

} // namespace Physics
