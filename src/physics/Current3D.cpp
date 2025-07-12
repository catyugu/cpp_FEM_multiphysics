#include "physics/Current3D.hpp"
#include <core/mesh/TetElement.hpp>
#include "utils/SimpleLogger.hpp"
#include "utils/Quadrature.hpp"
#include "utils/Exceptions.hpp"

namespace Physics {

Current3D::Current3D(const Core::Material& material) : material_(material) {}

const char* Current3D::getName() const { return "Electromagnetics 3D"; }
const char* Current3D::getVariableName() const { return "Voltage"; }

void Current3D::setup(Core::Mesh& mesh, Core::DOFManager& dof_manager) {
    mesh_ = &mesh;
    dof_manager_ = &dof_manager;
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Setting up ", getName(), " for mesh with material '", material_.getName(), "'.");

    size_t num_eq = dof_manager_->getNumEquations();
    K_.resize(num_eq, num_eq);
    F_.resize(num_eq, 1);
    U_.resize(num_eq, 1);
    U_prev_.resize(num_eq, 1);
    F_.setZero();
}

void Current3D::assemble() {
    auto& logger = SimpleLogger::Logger::instance();
    logger.info("Assembling system for ", getName());

    K_.setZero();
    // F_ is handled by applySources()

    double sigma = material_.getProperty("electrical_conductivity");
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity() * sigma;

    std::vector<Eigen::Triplet<double>> triplet_list;
    auto quadrature_points = Utils::Quadrature::getTetrahedronQuadrature(element_order_);

    int valid_elements_found = 0;
    for (const auto& elem_ptr : mesh_->getElements()) {
        auto* tet_elem = dynamic_cast<Core::TetElement*>(elem_ptr);
        if (tet_elem) {
            valid_elements_found++;
            Eigen::Matrix4d ke_local = Eigen::Matrix4d::Zero();
            double detJ = tet_elem->getVolume() * 6.0;
            for(const auto& qp : quadrature_points) {
                auto B = tet_elem->getBMatrix();
                ke_local += B.transpose() * D * B * qp.weight * detJ;
            }

            auto nodes = tet_elem->getNodes();
            int dofs[4];
            for(int i=0; i<4; ++i) {
                dofs[i] = dof_manager_->getEquationIndex(nodes[i]->getId(), getVariableName());
            }

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if (dofs[i] != -1 && dofs[j] != -1) {
                        triplet_list.emplace_back(dofs[i], dofs[j], ke_local(i, j));
                    }
                }
            }
        }
    }

    if (valid_elements_found == 0) {
        throw Exception::ConfigurationException(
            "Assembly failed for " + std::string(getName()) +
            ": No valid elements (TetElement) were found in the mesh. Check if the mesh is a 3D volume mesh."
        );
    }

    K_.setFromTriplets(triplet_list.begin(), triplet_list.end());

    logger.info("Assembly for ", getName(), " complete.");
}

} // namespace Physics
