#include "atmosphere.hpp"

Atmosphere::Atmosphere(/* args */)
{
    // Ficheiro com windspeed
}

double Atmosphere::getDensity(const double& altitude)
{
    return getPressure(altitude)/(R*getTemperature(altitude));
}

double Atmosphere::getPressure(const double& altitude)
{
    return pressure_SL*pow((1.0-0.0065*(altitude/temperature_SL)),5.2561);
}

double Atmosphere::getTemperature(const double& altitude)
{
    return temperature_SL - 6.5*(altitude/1000.0);
}

double Atmosphere::getSpeedOfSound(const double& altitude)
{
    return sqrt(gamma*R*getTemperature(altitude));
}

double Atmosphere::getDynamicViscosity(const double& altitude)
{
    return dViscosity_SL*((temperature_SL+113.0)/(getTemperature(altitude)+113.0))*(getTemperature(altitude)/temperature_SL);
}