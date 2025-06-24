#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Element.hpp"
#include "core/Node.hpp"
#include <vector>

namespace Physics {

Heat1D::Heat1D(const Core::Material& material) : material_(material), k_(0.0) {}

const char* Heat1D::getName() const { return "Heat Transfer 1D"; }
const char* Heat1D::getVariableName() const { return "Temperature"; }


void Heat1D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;

    k_ = material_.getProperty("thermal_conductivity");

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> Thermal Conductivity (k): ", k_);

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq);
    U_.resize(num_eq);

    volumetric_heat_source_.resize(mesh_->getElements().size(), 0.0);
}

void Heat1D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    std::vector<Eigen::Triplet<double>> triplet_list;

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        const auto& elem_ptr = mesh_->getElements()[i];
        Core::LineElement* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
        if (line_elem) {
            double h = line_elem->getLength();
            double ke_val = k_ / h;

            auto nodes = line_elem->getNodes();
            int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
            int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());

            triplet_list.emplace_back(dof_i, dof_i, ke_val);
            triplet_list.emplace_back(dof_i, dof_j, -ke_val);
            triplet_list.emplace_back(dof_j, dof_i, -ke_val);
            triplet_list.emplace_back(dof_j, dof_j, ke_val);

            double fe_val = volumetric_heat_source_[i] * h / 2.0;
            F_(dof_i) += fe_val;
            F_(dof_j) += fe_val;
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

void Heat1D::setVolumetricHeatSource(const std::vector<double>& source) {
    volumetric_heat_source_ = source;
    SimpleLogger::Logger::instance().info(getName(), ": Received volumetric heat source data.");
}

} // namespace Physics
