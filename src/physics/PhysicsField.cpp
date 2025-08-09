#include "physics/PhysicsField.hpp"
#include <core/mesh/LineElement.hpp>
#include <core/mesh/TetElement.hpp>
#include <core/mesh/TriElement.hpp>
#include <core/material/VariableManager.hpp>
#include "utils/SimpleLogger.hpp"
#include <set>

namespace Physics {
    void PhysicsField::setup(Core::Problem& problem, Core::Mesh& mesh, Core::DOFManager& dof_manager) {
        problem_ = &problem; // Store the pointer to the problem
        mesh_ = &mesh;
        dof_manager_ = &dof_manager;

        auto num_eq = static_cast<Eigen::Index>(dof_manager.getNumEquations());
        K_.resize(num_eq, num_eq);
        M_.resize(num_eq, num_eq);
        F_.resize(num_eq);
        F_coupling_.resize(num_eq); // 初始化耦合向量
        U_.resize(num_eq);
        U_prev_.resize(num_eq);

        F_.setZero();
        F_coupling_.setZero();
        U_.setZero();
        U_prev_.setZero();
    }

    // ... (其他函数保持不变，特别是 applyBCs) ...
    void PhysicsField::addBC(std::unique_ptr<Core::BoundaryCondition> bc) {
        bcs_.push_back(std::move(bc));
    }
    void PhysicsField::addBCs(std::vector<std::unique_ptr<Core::BoundaryCondition>> &&bcs) {
        for (auto& bc : bcs) {
            bcs_.push_back(std::move(bc));
        }
    }
    
    void PhysicsField::applyBCs() {
        auto &logger = Utils::Logger::instance();
        logger.info("Applying ", bcs_.size(), " defined BCs for ", getName());
        std::map<int, double> dirichlet_dofs;
        std::vector<const Core::BoundaryCondition *> other_bcs;
        for (const auto &bc: bcs_) {
            if (auto *dirichlet = dynamic_cast<const Core::DirichletBC *>(bc.get())) {
                if (dirichlet->getEquationIndex() != -1) {
                    dirichlet_dofs[dirichlet->getEquationIndex()] = dirichlet->getValue()(0);
                }
            } else {
                other_bcs.push_back(bc.get());
            }
        }
        logger.info("Applying ", dirichlet_dofs.size(), " unique, consolidated Dirichlet BCs.");
        if (!dirichlet_dofs.empty()) {
            Eigen::VectorXd U_bc = Eigen::VectorXd::Zero(K_.rows());
            for (const auto &pair: dirichlet_dofs) U_bc(pair.first) = pair.second;
            F_ -= K_ * U_bc;

            std::set<int> dirichlet_indices;
            for (const auto &pair: dirichlet_dofs) dirichlet_indices.insert(pair.first);

            for (int k = 0; k < K_.outerSize(); ++k) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(K_, k); it; ++it) {
                    if (dirichlet_indices.count(static_cast<int>(it.row())) || dirichlet_indices.count(static_cast<int>(it.col()))) {
                        it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
                    }
                }
            }
            K_.prune(0.0);
            for (const auto &pair: dirichlet_dofs) F_(pair.first) = pair.second;
        }
        for (const auto &bc: other_bcs) {
            bc->apply(K_, F_);
        }
    }

    void PhysicsField::applySources() {
        F_.setZero();
        for (const auto &source: source_terms_) {
            source->apply(F_, *dof_manager_, *mesh_, getVariableName(), element_order_);
        }
    }

    // ... (所有其他函数保持不变)
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
            Utils::Logger::instance().error("Cannot set initial conditions before field setup.");
        }
    }

    template<typename Func>
    void PhysicsField::setInitialConditions(Func func) {
        if (U_.size() > 0) {
            int num_components = getNumComponents();
            for (const auto &node: mesh_->getNodes()) {
                int base_dof_idx = dof_manager_->getEquationIndex(node->getId(), getVariableName());
                if (base_dof_idx != -1) {
                    const auto &coords = node->getCoords();
                    auto values = func(coords);
                    if constexpr (std::is_convertible_v<decltype(values), double>) {
                        // 修复：对于标量值，直接赋值给第一个分量
                        U_(base_dof_idx) = values;
                    } else {
                        if (values.size() == num_components) {
                            for (int c = 0; c < num_components; ++c) U_(base_dof_idx + c) = values(c);
                        }
                    }
                }
            }
            U_prev_ = U_;
        }
    }

    std::vector<int> PhysicsField::getElementDofs(Core::Element *elem) const {
        const auto &vertex_nodes = elem->getNodes();
        const size_t num_vertices = vertex_nodes.size();
        elem->setOrder(element_order_);
        const size_t num_elem_nodes = elem->getNumNodes();
        std::vector<int> base_dofs;
        base_dofs.reserve(num_elem_nodes);
        for (size_t i = 0; i < num_vertices; ++i) {
            base_dofs.push_back(dof_manager_->getEquationIndex(vertex_nodes[i]->getId(), getVariableName()));
        }

        if (element_order_ > 1) {
            if (dynamic_cast<const Core::LineElement *>(elem)) {
                base_dofs.push_back(
                    dof_manager_->getEdgeEquationIndex({vertex_nodes[0]->getId(), vertex_nodes[1]->getId()},
                                                       getVariableName()));
            } else if (dynamic_cast<const Core::TriElement *>(elem)) {
                const std::vector<std::pair<int, int> > edges = {{0, 1}, {1, 2}, {2, 0}};
                for (const auto &edge: edges) {
                    base_dofs.push_back(dof_manager_->getEdgeEquationIndex(
                        {vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName()));
                }
            } else if (dynamic_cast<const Core::TetElement *>(elem)) {
                const std::vector<std::pair<int, int> > edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                for (const auto &edge: edges) {
                    base_dofs.push_back(dof_manager_->getEdgeEquationIndex(
                        {vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName()));
                }
            }
        }

        std::vector<int> dofs;
        int num_components = getNumComponents();
        dofs.reserve(num_elem_nodes * num_components);
        for (int base_dof: base_dofs) {
            if (base_dof != -1) {
                for (int c = 0; c < num_components; ++c) dofs.push_back(base_dof + c);
            } else {
                for (int c = 0; c < num_components; ++c) dofs.push_back(-1);
            }
        }
        return dofs;
    }

    Eigen::SparseMatrix<double> &PhysicsField::getStiffnessMatrix() { return K_; }
    Eigen::SparseMatrix<double> &PhysicsField::getMassMatrix() { return M_; }
    Eigen::VectorXd &PhysicsField::getRHS() { return F_; }
    Eigen::VectorXd &PhysicsField::getCouplingRHS() { return F_coupling_; }
    Eigen::VectorXd &PhysicsField::getSolution() { return U_; }
    const Eigen::SparseMatrix<double> &PhysicsField::getStiffnessMatrix() const { return K_; }
    const Eigen::SparseMatrix<double> &PhysicsField::getMassMatrix() const { return M_; }
    const Eigen::VectorXd &PhysicsField::getRHS() const { return F_; }
    const Eigen::VectorXd &PhysicsField::getSolution() const { return U_; }

    // ========= 新增：变量管理和更新功能实现 =========

    void PhysicsField::updateElementVariables() {
        if (!mesh_ || !dof_manager_ || U_.size() == 0) {
            Utils::Logger::instance().warn("Cannot update element variables: missing mesh, DOF manager, or solution");
            return;
        }

        const std::string& var_name = getFieldVariableName();
        auto& logger = Utils::Logger::instance();
        logger.info("Updating element variables for field: ", getName(), " (variable: ", var_name, ")");

        // 遍历所有元素，更新每个元素的变量值
        for (auto& element : mesh_->getElements()) {
            // 获取这个元素的DOF索引
            std::vector<int> dof_indices = getElementDofs(element);

            // 从解向量中提取对应的值
            std::vector<double> solution_values;
            solution_values.reserve(dof_indices.size());

            for (int dof_idx : dof_indices) {
                if (dof_idx >= 0 && dof_idx < U_.size()) {
                    solution_values.push_back(U_(dof_idx));
                }
            }

            // 更新元素的变量值
            element->updateVariableFromSolution(var_name, solution_values, dof_indices);
        }

        logger.info("Element variables updated for ", mesh_->getElements().size(), " elements");
    }

    void PhysicsField::registerFieldVariable() {
        const std::string& var_name = getFieldVariableName();
        auto& var_manager = Core::VariableManager::getInstance();

        if (!var_manager.hasVariable(var_name)) {
            std::string description = std::string("Variable for ") + getName() + " physics field";
            var_manager.registerVariable(var_name, 0.0, description);

            auto& logger = Utils::Logger::instance();
            logger.info("Registered variable '", var_name, "' for physics field: ", getName());
        }
    }
}
