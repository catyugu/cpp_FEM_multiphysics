#ifndef SOLVER_H
#define SOLVER_H

#include "mesh.h"
#include <Eigen/Dense>

// Forward-declare HeatParams so we don't need to include HeatField.h here
struct HeatParams;

// --- Base Element Class ---
// An abstract interface for any type of finite element.
class FiniteElement {
public:
    virtual ~FiniteElement() = default;

    // A pure virtual function that each physics element must implement.
    virtual void calculate_local_matrices(const Node* const* element_nodes,
                                          const HeatParams& params,
                                          Eigen::MatrixXd& local_stiffness,
                                          Eigen::MatrixXd& local_mass) const = 0;
};

// --- Specific Physics Implementation ---
// Contains the specific formulas for a 4-node heat transfer tetrahedron.
class HeatTetrahedron : public FiniteElement {
public:
    void calculate_local_matrices(const Node* const* element_nodes,
                                  const HeatParams& params,
                                  Eigen::MatrixXd& local_stiffness,
                                  Eigen::MatrixXd& local_mass) const override;
};

#endif //SOLVER_H