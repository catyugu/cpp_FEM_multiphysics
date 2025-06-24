#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <string>
#include <vector>

// An enum to define the type of boundary condition
enum class BCType {
    Dirichlet, // Fixed value (e.g., temperature = 100 C)
    Neumann,   // Fixed flux (e.g., heat flux = 50 W/m^2)
    Robin      // Convection (e.g., h(T_s - T_inf))
};

// A struct to hold all data for a single BC
struct BoundaryCondition {
    // Identifies which part of the mesh this BC applies to
    std::vector<int> entity_ids; // Can be node IDs or face IDs
    std::string entity_type;       // "node" or "face"

    // Identifies which physics field this BC applies to
    std::string dof_name;          // e.g., "temperature", "voltage"

    // Defines the BC type and values
    BCType type;
    double value1{0.0}; // Used for Dirichlet value, Neumann flux, or Robin T_inf
    double value2{0.0}; // Used for Robin convection coefficient 'h'
};

#endif //BOUNDARY_CONDITIONS_H