#include "physics/Current2D.hpp"
#include <core/mesh/TriElement.hpp>
#include <core/mesh/Node.hpp>
#include "utils/SimpleLogger.hpp"
#include <vector>
#include "utils/Quadrature.hpp"

namespace Physics {

    Current2D::Current2D(const Core::Material& material)
        : material_(material) {}

    const char* Current2D::getName() const { return "Electromagnetics 2D"; }
    const char* Current2D::getVariableName() const { return "Voltage"; }

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
        auto quadrature_points = Utils::Quadrature::getTriangleQuadrature(element_order_);

        for (const auto& elem_ptr : mesh_->getElements()) {
            auto* tri_elem = dynamic_cast<Core::TriElement*>(elem_ptr);
            if (tri_elem) {
                double local_sigma = material_.getProperty("electrical_conductivity");
                Eigen::Matrix3d ke_local = Eigen::Matrix3d::Zero();
                double detJ = tri_elem->getArea() * 2.0;

                for(const auto& qp : quadrature_points) {
                    auto B = tri_elem->getBMatrix();
                    Eigen::Matrix2d D = Eigen::Matrix2d::Identity() * local_sigma;
                    ke_local += B.transpose() * D * B * qp.weight * detJ;
                }

                auto nodes = tri_elem->getNodes();
                int dofs[3];
                for(int i=0; i<3; ++i) dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
                for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) triplet_list.emplace_back(dofs[i], dofs[j], ke_local(i, j));
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
}