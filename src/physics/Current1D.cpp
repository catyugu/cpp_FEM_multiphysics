#include "physics/Current1D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Element.hpp"
#include "core/Node.hpp"

namespace Physics {

// --- FIX: Initialize heat_field_ to nullptr in the constructor ---
Current1D::Current1D(const Core::Material& material)
    : material_(material), heat_field_(nullptr) {}

const char* Current1D::getName() const { return "Electromagnetics 1D"; }
const char* Current1D::getVariableName() const { return "Voltage"; }
void Current1D::setCoupledHeatField(const PhysicsField* heat_field) { heat_field_ = heat_field; }

void Current1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    M_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    F_.setZero();
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);
}

void Current1D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    std::vector<Eigen::Triplet<double>> k_triplets;

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
        if (line_elem) {
            double T_avg = 300.0;
            if (heat_field_ && heat_field_->getSolution().size() > 0) {
                T_avg = 0.0;
                const auto& T_solution = heat_field_->getSolution();
                auto nodes = line_elem->getNodes();
                for(int i=0; i<2; ++i) {
                    int dof_idx = dof_manager_->getEquationIndex(nodes[i]->getId(), "Temperature");
                    T_avg += T_solution(dof_idx);
                }
                T_avg /= 2.0;
            }
            double local_sigma = material_.getProperty("electrical_conductivity", T_avg);
            double h = line_elem->getLength();
            double ke_val = local_sigma / h;

            auto nodes = line_elem->getNodes();
            int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
            int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());

            k_triplets.emplace_back(dof_i, dof_i, ke_val);
            k_triplets.emplace_back(dof_i, dof_j, -ke_val);
            k_triplets.emplace_back(dof_j, dof_i, -ke_val);
            k_triplets.emplace_back(dof_j, dof_j, ke_val);
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

std::vector<double> Current1D::calculateJouleHeat() const {
    std::vector<double> joule_heat(mesh_->getElements().size(), 0.0);
    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* line_elem = dynamic_cast<Core::LineElement*>(mesh_->getElements()[i]);
        if (line_elem) {
            double T_avg = 300.0;
            if (heat_field_ && heat_field_->getSolution().size() > 0) {
                 auto nodes = line_elem->getNodes();
                 T_avg = 0.5 * (heat_field_->getSolution()(dof_manager_->getEquationIndex(nodes[0]->getId(), "Temperature")) +
                                heat_field_->getSolution()(dof_manager_->getEquationIndex(nodes[1]->getId(), "Temperature")));
            }
            double local_sigma = material_.getProperty("electrical_conductivity", T_avg);
            double h = line_elem->getLength();
            auto nodes = line_elem->getNodes();
            int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
            int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());
            double V_i = U_(dof_i);
            double V_j = U_(dof_j);
            double E = std::abs(V_i - V_j) / h;
            joule_heat[i] = local_sigma * E * E;
        }
    }
    return joule_heat;
}
} // namespace Physics
