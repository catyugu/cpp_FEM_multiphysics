#ifndef PARAMS_H
#define PARAMS_H

#include "SimpleLogger.hpp"
struct Params {
    double k{300.0};   // Thermal conductivity of copper
    double rou{8940.0}; // Density of copper
    double Cp{385.0};  // Specific heat capacity of copper
};

#endif //PARAMS_H