#ifndef ROCKETSIM_HPP
#define ROCKETSIM_HPP

#include "f_maths.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

class RocketSim
{
private:
    std::vector<double> time, mass, thrust, drag;
public:
    RocketSim(/* args */);
    void timeMassThrustDragInput(const std::string& filename);
};

#endif
