#ifndef NODE_HPP
#define NODE_HPP

#include <vector>

namespace Core {

    class Node {
    public:
        // Constructor for a node with a given ID and coordinates
        Node(int id, double x, double y = 0.0, double z = 0.0);

        // Get the node ID
        int getId() const;

        // Get the node coordinates
        const std::vector<double>& getCoords() const;

    private:
        int id_;
        std::vector<double> coords_; // Store coordinates (x, y, z)
    };

} // namespace Core

#endif // NODE_HPP
