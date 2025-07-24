#ifndef PRESCRIBEDCURRENTDENSITY_HPP
#define PRESCRIBEDCURRENTDENSITY_HPP

#include "SourceTerm.hpp"
#include <Eigen/Dense>

namespace Core {

    class PrescribedCurrentDensity : public SourceTerm {
    public:
        PrescribedCurrentDensity(int element_id, const Eigen::Vector3d& current_density, const std::string& tag = "");

        void apply(Eigen::VectorXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const override;

    private:
        int element_id_;
        Eigen::Vector3d current_density_;
    };

} // namespace Core

#endif // PRESCRIBEDCURRENTDENSITY_HPP