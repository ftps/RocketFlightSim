#ifndef ROCKET_HPP
#define ROCKET_HPP

#include "f_maths.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

enum Comp{cone = 0, body = 1, motor = 2};
struct Component {
    double L, R, M, Xcm;
    uint type;

    Component(const std::string& aux_s, const double& Xtip, const uint& type);
    fm::Matrix33 momInertia();
};

class Rocket
{
public: 

    // State Variables
    fm::Vector3D X,L,P; //Position, Angular and Linear momentum
    fm::Quat Q; //Quartenion

    // Geometry Variables
    double Xrb, Rrb, Arb, ncant, rf, Xf, Xcm, Xcm_body, Xcm_motor;
    double lcc2, lcn2;
    uint nfin;
    std::vector<Component> comp;

    // Physics Variables
    fm::Matrix33 iTensorBody, iTensorMotor;
    double massBody, massMotor;
    std::vector<double> time, mass, thrust, drag, dmdt;

    // Coefficients
    double Cda, CA = 0, CN = 0, CR = 2*M_PI; // CR is for a flat plane

public:
    Rocket(const std::string& filename);
    void timeMassThrustDragInput(const std::string& filename);

    // Inercia calculation methods. Cada componente assume ser um cilindro com peso.
    fm::Matrix33 noseConeInercia(const double& noseConeLength,const double& noseConeRadius, const double& noseConeMass);
    fm::Matrix33 componentInercia(const double& componentLength,const double& componentRadius, const double& componentMass); 

    Rocket();

};

std::vector<std::string> spltString(const std::string& s, const char& cc = ' ');

#endif