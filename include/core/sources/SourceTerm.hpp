#ifndef SOURCETERM_HPP
#define SOURCETERM_HPP

#include <Eigen/Dense>
#include <string>
#include "core/DOFManager.hpp"
#include "core/mesh/Mesh.hpp"


namespace Core {

    class SourceTerm {
    public:
        explicit SourceTerm(std::string tag) : tag_(std::move(tag)) {}
        virtual ~SourceTerm() = default;

        // 这个函数签名现在对于编译器是完全清晰的
        // MODIFICATION: Added 'int element_order'
        virtual void apply(Eigen::VectorXd& F, const DOFManager& dof_manager, const Mesh& mesh, const std::string& var_name, int element_order) const = 0;

        const std::string& getTag() const { return tag_; }

    protected:
        std::string tag_;
    };

} // namespace Core

#endif // SOURCETERM_HPP