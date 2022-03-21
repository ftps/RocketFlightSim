#include "gnuplot-iostream.hpp"
#include "rocket.hpp"

double F(const double& t, const double& x)
{
    return -3*x;
}


int main(int argc, char* argv[])
{
    std::cout << "O Jogre é cringe" << std::endl;
    

    Gnuplot gp;
    
    fm::pair_tV<double> y, y2, y3;
    fm::function_tV<double> f = F;
    
    y = fm::ode::ODEsolveCT(f,0.0,1.0);

    for(double t = 0; t <= 10; t += 0.01){
        y2.emplace_back(t, exp(-3*t));
    }

    gp << "plot '-' with lines title 'cringe', '-' with lines title 'cringe2'\n";
    gp.send(y);
    gp.send(y2);


    return 0;
}