#include "physics/Current2D.hpp"
#include <core/mesh/TriElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/SimpleLogger.hpp"
#include <vector>

namespace Physics {

// Initialize heat_field_ to nullptr to prevent garbage pointer access
Current2D::Current2D(const Core::Material& material)
    : material_(material), heat_field_(nullptr) {}

const char* Current2D::getName() const { return "Electromagnetics 2D"; }
const char* Current2D::getVariableName() const { return "Voltage"; }
void Current2D::setCoupledHeatField(const PhysicsField* heat_field) { heat_field_ = heat_field; }

void Current2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
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

void Current2D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    F_.setZero();

    std::vector<Eigen::Triplet<double>> triplet_list;

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
        if (tri_elem) {

            double T_avg = 300.0;
            if (heat_field_ && heat_field_->getSolution().size() > 0) {
                const auto& T_solution = heat_field_->getSolution();
                auto nodes = tri_elem->getNodes();
                if (dof_manager_->getEquationIndex(nodes[0]->getId(), "Temperature") == -1 ||
                    dof_manager_->getEquationIndex(nodes[1]->getId(), "Temperature") == -1 ||
                    dof_manager_->getEquationIndex(nodes[2]->getId(), "Temperature") == -1) {
                    logger.warn("    Element ", tri_elem->getId(), " has missing temperature DOFs. Using default T_avg.");
                } else {
                    T_avg = (T_solution(dof_manager_->getEquationIndex(nodes[0]->getId(), "Temperature")) +
                             T_solution(dof_manager_->getEquationIndex(nodes[1]->getId(), "Temperature")) +
                             T_solution(dof_manager_->getEquationIndex(nodes[2]->getId(), "Temperature"))) / 3.0;
                }
            }

            double local_sigma = material_.getProperty("electrical_conductivity", T_avg);
            double area = tri_elem->getArea();
            if (area <= 0) {
                logger.error("    Element ", tri_elem->getId(), " has zero or negative area! Skipping.");
                continue;
            }

            auto B = tri_elem->getBMatrix();
            Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * local_sigma;
            Eigen::Matrix3d ke = area * B.transpose() * D * B;

            auto nodes = tri_elem->getNodes();
            int dofs[3];
            for(int i=0; i<3; ++i) dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) triplet_list.emplace_back(dofs[i], dofs[j], ke(i, j));
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
std::vector<double> Current2D::calculateJouleHeat() const {
    std::vector<double> joule_heat(mesh_->getElements().size(), 0.0);
    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(mesh_->getElements()[i]);
        if (tri_elem) {
            double T_avg = 300.0;
            if (heat_field_ && heat_field_->getSolution().size() > 0) {
                T_avg = 0.0;
                const auto& T_solution = heat_field_->getSolution();
                auto nodes = tri_elem->getNodes();
                for(int j=0; j<3; ++j) {
                    int dof_idx = dof_manager_->getEquationIndex(nodes[j]->getId(), "Temperature");
                    T_avg += T_solution(dof_idx);
                }
                T_avg /= 3.0;
            }
            double local_sigma = material_.getProperty("electrical_conductivity", T_avg);
            auto B = tri_elem->getBMatrix();
            auto nodes = tri_elem->getNodes();
            Eigen::Vector3d nodal_voltages;
            for(int j=0; j<3; ++j) {
                int dof_idx = dof_manager_->getEquationIndex(nodes[j]->getId(), getVariableName());
                nodal_voltages(j) = U_(dof_idx);
            }
            Eigen::Vector2d grad_V = B * nodal_voltages;
            joule_heat[i] = local_sigma * grad_V.squaredNorm();
        }
    }
    return joule_heat;
}

} // namespace Physics
