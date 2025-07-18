#ifndef VOLUMETRICSOURCE_HPP
#define VOLUMETRICSOURCE_HPP

#include "SourceTerm.hpp"
#include <string>

namespace Core {
    class Mesh;
    class DOFManager;

    class VolumetricSource : public SourceTerm {
    public:
        // Constructor is now simpler
        VolumetricSource(int element_id, double total_power, const std::string& tag);
        void apply(Eigen::MatrixXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const override;

    private:
        int element_id_;
        double total_power_;
    };

} // namespace Core

#endif // VOLUMETRICSOURCE_HPP