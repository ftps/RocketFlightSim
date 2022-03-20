#ifndef ROCKET_HPP
#define ROCKET_HPP

#include <vector>

class Rocket
{
public: 

    // State Variables
    std::vector<double> X,Q,L,P; //Position, Quartenion, Angular and Linear momentum


    Rocket();

};



#endif