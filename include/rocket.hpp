#ifndef ROCKET_HPP
#define ROCKET_HPP

#include "f_maths.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

class Rocket
{
public: 

    // State Variables
    fm::Vector3D X,L,P; //Position, Angular and Linear momentum
    fm::Quat Q; //Quartenion

    // Geometry Variables
    double noseConeLength, noseConeRadius, noseConeMass;
    double tubeLength, tubeRadius;

    // Physics Variables
    fm::Matrix33 inerciaTensor;
    std::vector<double> time, mass, thrust, drag;

public:

    void timeMassThrustDragInput(const std::string& filename);

    // Inercia calculation methods. Cada componente assume ser um cilindro com peso.
    void noseConeInercia(const double& noseConeLength,const double& noseConeRadius, const double& noseConeMass);
    void componentInercia(const double& componentLength,const double& componentRadius, const double& componentMass); 

    Rocket();

};



#endif