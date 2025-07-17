#ifndef ELEMENTGEOMETRY_HPP
#define ELEMENTGEOMETRY_HPP

#include <vector>
#include <memory>
#include <Eigen/Dense>
#include "core/mesh/Node.hpp"

namespace Core {

    /**
     * @class ElementGeometry
     * @brief Holds the raw geometric information of an element's vertices.
     */
    class ElementGeometry {
    public:
        // Constructor now takes the dimension explicitly.
        explicit ElementGeometry(const std::vector<Node*>& vertex_nodes, int dimension);

        const Eigen::MatrixXd& get_vertex_coords() const { return vertex_coords_; }
        int get_dimension() const { return dimension_; }
        size_t get_num_vertices() const { return vertex_nodes_.size(); }

    private:
        std::vector<Node*> vertex_nodes_;
        Eigen::MatrixXd vertex_coords_;
        int dimension_;
    };

} // namespace Core

#endif // ELEMENTGEOMETRY_HPP