#include "physics/EMag2D.hpp"
#include "core/TriElement.hpp"
#include "core/Node.hpp"
#include "utils/SimpleLogger.hpp"

namespace Physics {

EMag2D::EMag2D(const Core::Material& material) : material_(material) {}

const char* EMag2D::getName() const { return "Electromagnetics 2D"; }
const char* EMag2D::getVariableName() const { return "Voltage"; }

void EMag2D::setCoupledHeatField(const PhysicsField* heat_field) {
    heat_field_ = heat_field;
}

void EMag2D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;

    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq);
    F_.setZero();
    U_.resize(num_eq);
    U_.setZero();
}

void EMag2D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    if (!heat_field_) {
        logger.warn(getName(), " is not coupled to a heat field. Assuming constant properties.");
    }

    K_.setZero();
    F_.setZero();

    std::vector<Eigen::Triplet<double>> triplet_list;

    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
        if (tri_elem) {
            // --- Temperature-dependent conductivity calculation ---
            double T_avg = 300.0; // Default temperature if uncoupled
            if (heat_field_) {
                T_avg = 0.0;
                const auto& T_solution = heat_field_->getSolution();
                auto nodes = tri_elem->getNodes();
                for(int i=0; i<3; ++i) {
                    int dof_idx = dof_manager_->getEquationIndex(nodes[i]->getId(), "Temperature");
                    T_avg += T_solution(dof_idx);
                }
                T_avg /= 3.0;
            }

            double local_sigma = material_.getProperty("electrical_conductivity", T_avg);
            Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * local_sigma;

            // --- Assembly using local conductivity ---
            double area = tri_elem->getArea();
            auto B = tri_elem->getBMatrix();
            Eigen::Matrix3d ke = area * B.transpose() * D * B;

            auto nodes = tri_elem->getNodes();
            int dofs[3];
            for(int i=0; i<3; ++i) {
                dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            }

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    triplet_list.emplace_back(dofs[i], dofs[j], ke(i, j));
                }
            }
        }
    }
    K_.setFromTriplets(triplet_list.begin(), triplet_list.end());

    // Stabilize matrix for other physics
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

std::vector<double> EMag2D::calculateJouleHeat() const {
    std::vector<double> joule_heat(mesh_->getElements().size(), 0.0);

    for (size_t i = 0; i < mesh_->getElements().size(); ++i) {
        auto* tri_elem = dynamic_cast<Core::TriElement*>(mesh_->getElements()[i]);
        if (tri_elem) {
            // --- Get local conductivity for heat calculation ---
             double T_avg = 300.0;
            if (heat_field_) {
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

            // --- Calculate Joule Heat ---
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
