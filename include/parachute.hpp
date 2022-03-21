#ifndef PARACHUTE_HPP
#define PARACHUTE_HPP

#include "f_maths.hpp"

class Parachute
{
public: 

    // State Variables
    fm::Vector3D X,P; //Position and Linear momentum

    // Parachute Parameters
    double chute_drag, chute_area; // Parachute Drag Coeff and area.

    Parachute();

};



#endif