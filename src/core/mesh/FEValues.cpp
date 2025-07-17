#include "core/FEValues.hpp"
#include "utils/ShapeFunctions.hpp"
#include "utils/Quadrature.hpp"
#include "core/mesh/Node.hpp"
#include "core/mesh/LineElement.hpp"  // Include derived element types
#include "core/mesh/TriElement.hpp"
#include "core/mesh/TetElement.hpp"
#include <stdexcept>

namespace Core {

// FIX: Constructor now takes a const Element&
FEValues::FEValues(const Element& elem, int order, int quad_order)
    : element_(elem), element_order_(order), quadrature_order_(quad_order) {

    // Get geometry and nodes from the element
    const ElementGeometry& geom = element_.getGeometry();
    const auto& nodes = geom.getNodes();
    const int dim = nodes[0]->getCoords().size();

    // This calculation for num_elem_nodes is simplified and may need adjustment
    // for a fully generic implementation, but works for this structure.
    size_t num_elem_nodes;
    if (dynamic_cast<const LineElement*>(&element_)) {
        num_elem_nodes = order + 1;
    } else if (dynamic_cast<const TriElement*>(&element_)) {
        num_elem_nodes = (order + 1) * (order + 2) / 2;
    } else if (dynamic_cast<const TetElement*>(&element_)) {
        num_elem_nodes = (order + 1) * (order + 2) * (order + 3) / 6;
    } else {
        throw std::runtime_error("Unsupported element type in FEValues");
    }


    // FIX: Determine Quadrature rule based on the element's actual type
    if (dynamic_cast<const LineElement*>(&element_)) {
        quadrature_points_ = Utils::Quadrature::getLineQuadrature(quad_order);
    } else if (dynamic_cast<const TriElement*>(&element_)) {
        quadrature_points_ = Utils::Quadrature::getTriangleQuadrature(quad_order);
    } else if (dynamic_cast<const TetElement*>(&element_)) {
        quadrature_points_ = Utils::Quadrature::getTetrahedronQuadrature(quad_order);
    } else {
        throw std::runtime_error("Unsupported element type in FEValues for Quadrature");
    }

    // Pre-calculate and cache values for all quadrature points
    for (const auto& qp : quadrature_points_) {
        // 1. Get shape functions and their derivatives in natural coordinates
        Eigen::VectorXd N_nat;
        Eigen::MatrixXd dN_dnat;

        // This logic needs to be robust for different element types
        // Example for Triangle, you would add Line and Tet cases
        if (dynamic_cast<const TriElement*>(&element_)) {
             N_nat = Utils::ShapeFunctions::getTriShapeFunctions(element_order_, qp.point(0), qp.point(1));
             dN_dnat = Utils::ShapeFunctions::getTriShapeFunctionDerivatives(element_order_, qp.point(0), qp.point(1));
        } else {
             // Add cases for LineElement and TetElement here
             throw std::runtime_error("Shape function logic not implemented for this element type in FEValues");
        }


        // 2. Calculate Jacobian matrix J = sum(node_coords * dN_dnat)
        Eigen::MatrixXd node_coords(dim, nodes.size());
        for(size_t i = 0; i < nodes.size(); ++i) {
            const auto& coords = nodes[i]->getCoords();
            for(int j=0; j<dim; ++j) {
                node_coords(j, i) = coords[j];
            }
        }

        Eigen::MatrixXd J = node_coords * dN_dnat;
        double detJ = J.determinant();

        // 3. Calculate gradients in real coordinates: dN/dx = dN/d_xi * J^-1
        Eigen::MatrixXd dNdx = dN_dnat * J.inverse();


        // 4. Cache the results
        N_values_.push_back(N_nat);
        dNdx_values_.push_back(dNdx.transpose()); // B matrix is often the transpose of this
        detJ_x_w_values_.push_back(detJ * qp.weight);
    }
}

void FEValues::reinit(int q_point_index) {
    if (q_point_index < 0 || q_point_index >= (int)N_values_.size()) {
        throw std::out_of_range("Quadrature point index is out of bounds.");
    }
    N_ = N_values_[q_point_index];
    dNdx_ = dNdx_values_[q_point_index];
    detJ_x_w_ = detJ_x_w_values_[q_point_index];
}

const Eigen::VectorXd& FEValues::get_shape_values() const {
    return N_;
}

const Eigen::MatrixXd& FEValues::get_shape_gradients() const {
    return dNdx_;
}

double FEValues::get_detJ_times_weight() const {
    return detJ_x_w_;
}

} // namespace Core