#ifndef ROCKET_HPP
#define ROCKET_HPP

#include "f_maths.hpp"

class Rocket
{
public: 

    // State Variables
    fm::Vector3D X,L,P; //Position, Angular and Linear momentum
    fm::Quat Q; //Quartenion


    Rocket();

};



#endif