#include "physics/PhysicsField.hpp"

#include <core/mesh/LineElement.hpp>
#include <core/mesh/TetElement.hpp>
#include <core/mesh/TriElement.hpp>

#include "utils/SimpleLogger.hpp"

namespace Physics {
    void PhysicsField::addBC(std::unique_ptr<Core::BoundaryCondition> bc) {
        bcs_.push_back(std::move(bc));
    }

    void PhysicsField::applyBCs() {
        auto &logger = Utils::Logger::instance();
        logger.info("Applying ", bcs_.size(), " defined BCs for ", getName());

        std::map<int, const Core::BoundaryCondition *> consolidated_bcs;
        for (const auto &bc: bcs_) {
            if (bc->getEquationIndex() != -1) {
                consolidated_bcs[bc->getEquationIndex()] = bc.get();
            }
        }

        logger.info("Applying ", consolidated_bcs.size(), " unique, consolidated BCs.");

        // Iterate through the consolidated map and apply the final, unique BCs.
        for (const auto &pair: consolidated_bcs) {
            pair.second->apply(K_, F_);
        }
    }


    void PhysicsField::applySources() {
        // Clear the force vector before applying sources to avoid accumulation across steps/iterations
        F_.setZero();
        for (const auto &source: source_terms_) {
            source->apply(F_, *dof_manager_, *mesh_, getVariableName());
        }
    }


    void PhysicsField::removeBCsByTag(const std::string &tag) {
        if (tag.empty()) return; // Do not remove BCs with no tag
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


    const Eigen::MatrixXd &PhysicsField::getPreviousSolution() const {
        return U_prev_;
    }

    void PhysicsField::setInitialConditions(double initial_value) {
        if (U_.size() > 0) {
            U_.setConstant(initial_value);
            U_prev_ = U_;
            Utils::Logger::instance().info("Set initial condition for '", getVariableName(), "' to ", initial_value);
        } else {
            Utils::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(),
                                            "'.");
        }
    }

    template<typename F>
    void PhysicsField::setInitialConditions(std::function<F> func) {
        if (U_.size() > 0) {
            for (int i = 0; i < U_.size(); ++i) {
                U_(i) = func(mesh_->getNode(i));
            }
            U_prev_ = U_;
            Utils::Logger::instance().info("Set initial condition for '", getVariableName(), "' to ", func);
        } else {
            Utils::Logger::instance().error("Cannot set initial conditions before field setup for '", getVariableName(),
                                            "'.");
        }
    }

  std::vector<int> PhysicsField::get_element_dofs(Core::Element* elem) const {
        const auto& vertex_nodes = elem->getNodes();
        const size_t num_vertices = vertex_nodes.size();
        elem->setOrder(element_order_);
        const size_t num_elem_nodes = elem->getNumNodes();
        std::vector<int> dofs(num_elem_nodes);

        // --- 1. Get Vertex DOFs (always first) ---
        for (size_t i = 0; i < num_vertices; ++i) {
            dofs[i] = dof_manager_->getEquationIndex(vertex_nodes[i]->getId(), getVariableName());
        }

        // --- 2. Get Higher-Order DOFs (if any) in their correct canonical order ---
        if (element_order_ > 1) {
            if (auto* line = dynamic_cast<const Core::LineElement*>(elem)) {
                // Canonical order for P2 Line: [v0, midpoint, v1]
                dofs[2] = dofs[1]; // Temporarily move v1's DOF to the end
                dofs[1] = dof_manager_->getEdgeEquationIndex({vertex_nodes[0]->getId(), vertex_nodes[1]->getId()}, getVariableName());
            }
            else if (auto* tri = dynamic_cast<const Core::TriElement*>(elem)) {
                // Canonical edge order for Tri6: v0,v1,v2, edge(0,1), edge(1,2), edge(2,0)
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}};
                int edge_dof_idx = 3;
                for (const auto& edge : edges) {
                    dofs[edge_dof_idx++] = dof_manager_->getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName());
                }
            } else if (auto* tet = dynamic_cast<const Core::TetElement*>(elem)) {
                // Canonical edge order for Tet10: v0,v1,v2,v3, edge(0,1), edge(0,2), edge(0,3), edge(1,2), edge(1,3), edge(2,3)
                const std::vector<std::pair<int, int>> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
                int edge_dof_idx = 4;
                for (const auto& edge : edges) {
                    dofs[edge_dof_idx++] = dof_manager_->getEdgeEquationIndex({vertex_nodes[edge.first]->getId(), vertex_nodes[edge.second]->getId()}, getVariableName());
                }
            } else {
                throw std::runtime_error("get_element_dofs: Unsupported element type for higher orders.");
            }
        }
        return dofs;
    }
    // --- END OF FIX ---
    Eigen::SparseMatrix<double> &PhysicsField::getStiffnessMatrix() { return K_; }
    Eigen::SparseMatrix<double> &PhysicsField::getMassMatrix() { return M_; }
    Eigen::MatrixXd &PhysicsField::getRHS() { return F_; }
    Eigen::MatrixXd &PhysicsField::getSolution() { return U_; }

    const Eigen::SparseMatrix<double> &PhysicsField::getStiffnessMatrix() const { return K_; }
    const Eigen::SparseMatrix<double> &PhysicsField::getMassMatrix() const { return M_; }
    const Eigen::MatrixXd &PhysicsField::getRHS() const { return F_; }
    const Eigen::MatrixXd &PhysicsField::getSolution() const { return U_; }
} // namespace Physics