#ifndef PARACHUTE_HPP
#define PARACHUTE_HPP

#include <vector>

class Parachute
{
public: 

    // State Variables
    std::vector<double> X,P; //Position and Linear momentum

    // Parachute Parameters
    double chute_drag, chute_area; // Parachute Drag Coeff and area.

    Parachute();

};



#endif