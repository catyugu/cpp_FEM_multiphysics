#include "physics/Current1D.hpp"
#include "utils/SimpleLogger.hpp"
#include <core/mesh/Element.hpp>
#include <core/mesh/Node.hpp>
#include "utils/Quadrature.hpp"

namespace Physics {
    Current1D::Current1D(const Core::Material& material)
        : material_(material) {}

    const char* Current1D::getName() const { return "Electromagnetics 1D"; }
    const char* Current1D::getVariableName() const { return "Voltage"; }
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
        auto quadrature_points = Utils::Quadrature::getLineQuadrature(element_order_);

        for (const auto& elem_ptr : mesh_->getElements()) {
            auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
            if (line_elem) {
                double local_sigma = material_.getProperty("electrical_conductivity");
                double h = line_elem->getLength();

                Eigen::Matrix2d ke_local = Eigen::Matrix2d::Zero();
                for(const auto& qp : quadrature_points) {
                    // For linear elements, B is constant, and detJ is length/2
                    double detJ = h / 2.0;
                    Eigen::Matrix<double, 1, 2> B;
                    B << -1/h, 1/h;
                    ke_local += B.transpose() * local_sigma * B * qp.weight * detJ;
                }

                auto nodes = line_elem->getNodes();
                int dof_i = dof_manager_->getEquationIndex(nodes[0]->getId(), getVariableName());
                int dof_j = dof_manager_->getEquationIndex(nodes[1]->getId(), getVariableName());
                int dofs[2] = {dof_i, dof_j};

                for(int i = 0; i < 2; ++i) {
                    for(int j = 0; j < 2; ++j) {
                        k_triplets.emplace_back(dofs[i], dofs[j], ke_local(i,j));
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
}