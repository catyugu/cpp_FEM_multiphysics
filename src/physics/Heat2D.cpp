#include "physics/Heat2D.hpp"
#include "core/TriElement.hpp"
#include "core/Node.hpp"
#include "utils/SimpleLogger.hpp"

namespace Physics {

Heat2D::Heat2D(const Core::Material& material) : material_(material), k_(0.0) {}

const char* Heat2D::getName() const { return "Heat Transfer 2D"; }
const char* Heat2D::getVariableName() const { return "Temperature"; }

void Heat2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    k_ = material_.getProperty("thermal_conductivity");

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");
    logger.info("-> Thermal Conductivity (k): ", k_);

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq);
    F_.setZero();
    U_.resize(num_eq);
    U_.setZero();
    volumetric_heat_source_.resize(mesh_->getElements().size(), 0.0);
}

void Heat2D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * k_;
    std::vector<Eigen::Triplet<double>> triplet_list;

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        const auto& elem_ptr = mesh_->getElements()[i];
        auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
        if (tri_elem) {
            double area = tri_elem->getArea();
            auto B = tri_elem->getBMatrix();
            Eigen::Matrix3d ke = area * B.transpose() * D * B;

            auto nodes = tri_elem->getNodes();
            int dofs[3];
            for(int j=0; j<3; ++j) {
                dofs[j] = dof_manager_->getEquationIndex(nodes[j]->getId(), getVariableName());
            }

            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    triplet_list.emplace_back(dofs[r], dofs[c], ke(r, c));
                }
            }

            double fe_val = volumetric_heat_source_[i] * area / 3.0;
            if (fe_val != 0.0) {
                for (int j=0; j<3; ++j) {
                    F_(dofs[j]) += fe_val;
                }
            }
        }
    }
    K_.setFromTriplets(triplet_list.begin(), triplet_list.end());

    // --- FIX: Ensure non-active DOFs don't create a singular matrix ---
    const std::string my_var = getVariableName();
    for(const auto& var_name : dof_manager_->getVariableNames()) {
        if (var_name != my_var) {
            // This is another physics field's variable. Set its diagonal to 1.
            for(const auto& node : mesh_->getNodes()) {
                int dof_idx = dof_manager_->getEquationIndex(node->getId(), var_name);
                if (dof_idx != -1) {
                    K_.coeffRef(dof_idx, dof_idx) = 1.0;
                }
            }
        }
    }

    logger.info("Assembly for ", getName(), " complete.");
}

void Heat2D::setVolumetricHeatSource(const std::vector<double>& source) {
    if (source.size() != mesh_->getElements().size()) {
        SimpleLogger::Logger::instance().error("Heat source vector size mismatch in Heat2D.");
        return;
    }
    volumetric_heat_source_ = source;
    SimpleLogger::Logger::instance().info(getName(), ": Received volumetric heat source data.");
}

} // namespace Physics
