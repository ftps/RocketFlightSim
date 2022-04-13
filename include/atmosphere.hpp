#ifndef ATMOSPHERE_HPP
#define ATMOSPHERE_HPP

#include "f_maths.hpp"


// Válido apenas para baixo da Tropopausa, ISA standart
class Atmosphere
{
private:
    // In function of Z
    std::vector<double> windspeed;

    // Constants at sea level
    double pressure_SL = 101325;
    double density_SL = 1.225;
    double temperature_SL = 288.15;
    double speed_of_sound_SL = 350.294;
    double gravity_SL = 9.80665;
    double R = 8.314462;
    double gamma = 1.4;
    double dViscosity_SL = 0.0000789;

public:
    Atmosphere(/* args */);
    double getTemperature(const double& altitude);
    double getPressure(const double& altitude);
    double getDensity(const double& altitude);
    double getSpeedOfSound(const double& altitude);
    double getDynamicViscosity(const double& altitude); //Sutherland Model, C = 113

};



#endif