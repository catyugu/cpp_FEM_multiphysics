#ifndef HEAT_FIELD_H
#define HEAT_FIELD_H

#include "SimpleLogger.hpp"
struct HeatParams {
    double k{300.0};   // Thermal conductivity of copper
    double rou{8940.0}; // Density of copper
    double Cp{385.0};  // Specific heat capacity of copper
};

#endif //HEAT_FIELD_H