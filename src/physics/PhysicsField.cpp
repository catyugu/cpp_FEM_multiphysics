#include "physics/PhysicsField.hpp"
#include <core/mesh/LineElement.hpp>
#include <core/mesh/TetElement.hpp>
#include <core/mesh/TriElement.hpp>
#include "utils/SimpleLogger.hpp"
#include <set> // 确保包含 set

namespace Physics {
    void PhysicsField::addBC(std::unique_ptr<Core::BoundaryCondition> bc) {
        bcs_.push_back(std::move(bc));
    }

    void PhysicsField::addBCs(std::vector<std::unique_ptr<Core::BoundaryCondition>> bcs) {
        bcs_.insert(bcs_.end(), std::make_move_iterator(bcs.begin()), std::make_move_iterator(bcs.end()));
    }

    void PhysicsField::applyBCs() {
        auto &logger = Utils::Logger::instance();
        logger.info("Applying ", bcs_.size(), " defined BCs for ", getName());

        // --- 高效的批量边界条件处理 ---
        std::map<int, double> dirichlet_dofs;
        std::vector<const Core::BoundaryCondition*> other_bcs;

        // 1. 分离狄利克雷和其他边界条件
        for (const auto& bc : bcs_) {
            if (auto* dirichlet = dynamic_cast<const Core::DirichletBC*>(bc.get())) {
                if (dirichlet->getEquationIndex() != -1) {
                    // 使用 .at(0) 因为我们的BC值都是1维向量
                    dirichlet_dofs[dirichlet->getEquationIndex()] = dirichlet->getValue()(0);
                }
            } else {
                other_bcs.push_back(bc.get());
            }
        }

        logger.info("Applying ", dirichlet_dofs.size(), " unique, consolidated Dirichlet BCs.");

        // 2. 批量修改 F 向量
        Eigen::VectorXd U_bc = Eigen::VectorXd::Zero(K_.rows());
        for (const auto& pair : dirichlet_dofs) {
            U_bc(pair.first) = pair.second;
        }
        F_ -= K_ * U_bc;

        // 3. 批量修改 K 矩阵
        std::set<int> dirichlet_indices;
        for(const auto& pair : dirichlet_dofs) {
            dirichlet_indices.insert(pair.first);
        }

        for (int k=0; k < K_.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(K_, k); it; ++it) {
                if (dirichlet_indices.count(it.row()) || dirichlet_indices.count(it.col())) {
                    if (it.row() == it.col()) {
                        it.valueRef() = 1.0; // 对角线元素设为1
                    } else {
                        it.valueRef() = 0.0; // 非对角线元素设为0
                    }
                }
            }
        }
        K_.prune(0.0); // 移除被置为0的元素

        // 4. 最终设置 F 向量
        for (const auto& pair : dirichlet_dofs) {
            F_(pair.first) = pair.second;
        }

        // 5. 应用其他类型的边界条件 (Neumann, Cauchy)
        for (const auto& bc : other_bcs) {
            bc->apply(K_, F_);
        }
    }

    void PhysicsField::applySources() {
        F_.setZero();
        for (const auto &source: source_terms_) {
            source->apply(F_, *dof_manager_, *mesh_, getVariableName(), element_order_);
        }
    }

    void PhysicsField::removeBCsByTag(const std::string &tag) {
        if (tag.empty()) return;
        auto it = std::remove_if(bcs_.begin(), bcs_.end(),
                                 [&](const std::unique_ptr<Core::BoundaryCondition> &bc_ptr) {
                                     return bc_ptr->getTag() == tag;
                                 });
        if (it != bcs_.end()) {
            bcs_.erase(it, bcs_.end());
        }
    }

    void PhysicsField::addSource(std::unique_ptr<Core::SourceTerm> source) {
        source_terms_.push_back(std::move(source));
    }

    void PhysicsField::removeSourcesByTag(const std::string &tag) {
        auto it = std::remove_if(source_terms_.begin(), source_terms_.end(),
                                 [&](const std::unique_ptr<Core::SourceTerm> &source_ptr) {
                                     return source_ptr->getTag() == tag;
                                 });
        if (it != source_terms_.end()) {
            source_terms_.erase(it, source_terms_.end());
        }
    }

    const std::vector<std::unique_ptr<Core::BoundaryCondition> > &PhysicsField::getBCs() const {
        return bcs_;
    }

    void PhysicsField::updatePreviousSolution() {
        U_prev_ = U_;
    }

    void PhysicsField::setElementOrder(int order) {
        if (order < 1) {
            Utils::Logger::instance().error("Element order must be at least 1.");
            return;
        }
        element_order_ = order;
    }

    const Eigen::VectorXd &PhysicsField::getPreviousSolution() const {
        return U_prev_;
    }

    void PhysicsField::setInitialConditions(double initial_value) {
        if (U_.size() > 0) {
            U_.setConstant(initial_value);
            U_prev_ = U_;
        } else {
            Utils::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(), "'.");
        }
    }

    template<typename F>
    void PhysicsField::setInitialConditions(std::function<F> func) {
        if (U_.size() > 0) {
            for (int i = 0; i < U_.size(); ++i) {
                U_(i) = func(mesh_->getNode(i));
            }
            U_prev_ = U_;
        } else {
            Utils::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(), "'.");
        }
    }

    std::vector<int> PhysicsField::getElementDofs(Core::Element* elem) const {
        const auto& vertex_nodes = elem->getNodes();
        const size_t num_vertices = vertex_nodes.size();
        elem->setOrder(element_order_);
        const size_t num_elem_nodes = elem->getNumNodes();
        std::vector<int> base_dofs;
        base_dofs.reserve(num_elem_nodes);

        for (size_t i = 0; i < num_vertices; ++i) {
            base_dofs.push_back(dof_manager_->getEquationIndex(vertex_nodes[i]->getId(), getVariableName()));
        }

        if (element_order_ > 1) {
            if (dynamic_cast<const Core::LineElement*>(elem)) {
                base_dofs.push_back(dof_manager_->getEdgeEquationIndex({vertex_nodes[0]->getId(), vertex_nodes[1]->getId()}, getVariableName()));
            } else if (dynamic_cast<const Core::TriElement*>(elem)) {
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto& edge : edges) {
                    base_dofs.push_back(dof_manager_->getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName()));
                }
            } else if (dynamic_cast<const Core::TetElement*>(elem)) {
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                for (const auto& edge : edges) {
                    base_dofs.push_back(dof_manager_->getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName()));
                }
            }
        }

        std::vector<int> dofs;
        int num_components = getNumComponents();
        dofs.reserve(num_elem_nodes * num_components);
        for (int base_dof : base_dofs) {
            if (base_dof != -1) {
                for (int c = 0; c < num_components; ++c) {
                    dofs.push_back(base_dof + c);
                }
            } else {
                for (int c = 0; c < num_components; ++c) {
                    dofs.push_back(-1);
                }
            }
        }
        return dofs;
    }

    Eigen::SparseMatrix<double> &PhysicsField::getStiffnessMatrix() { return K_; }
    Eigen::SparseMatrix<double> &PhysicsField::getMassMatrix() { return M_; }
    Eigen::VectorXd &PhysicsField::getRHS() { return F_; }
    Eigen::VectorXd &PhysicsField::getSolution() { return U_; }

    const Eigen::SparseMatrix<double> &PhysicsField::getStiffnessMatrix() const { return K_; }
    const Eigen::SparseMatrix<double> &PhysicsField::getMassMatrix() const { return M_; }
    const Eigen::VectorXd &PhysicsField::getRHS() const { return F_; }
    const Eigen::VectorXd &PhysicsField::getSolution() const { return U_; }
}