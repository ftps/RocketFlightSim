#include "gnuplot-iostream.hpp"
#include "rocket.hpp"
#include "rocketSim.hpp"

double F(const double& t, const double& x)
{
    return -3*x;
}


int main(int argc, char* argv[])
{
    //std::cout << "O Jogre é cringe" << std::endl;
    

    // Gnuplot gp;
    
    // fm::pair_tV<double> y, y2, y3;
    // fm::function_tV<double> f = F;
    
    // y = fm::ode::ODEsolveCT(f,0.0,1.0);

    // for(double t = 0; t <= 10; t += 0.01){
    //     y2.emplace_back(t, exp(-3*t));
    // }

    // gp << "plot '-' with lines title 'cringe', '-' with lines title 'cringe2'\n";
    // gp.send(y);
    // gp.send(y2);

    // Original data
//    std::vector<double> xData = { 1, 5, 10, 15, 20 };
//    std::vector<double> yData = { 0.3, 0.5, 0.8, 0.1, 0.14 };

//    // Set up some points for interpolation in xVals
//    const int NPTS = 20;
//    std::vector<double> xVals, yVals;
//    for ( int i = 1; i <= NPTS; i++ ) xVals.push_back( (double)i );

//    // Interpolate
//    for ( double x : xVals )
//    {
//       double y = fm::interp::polyInterp( xData, yData, x );
//       yVals.push_back( y );
//    }

// // Output
//    #define SP << std::fixed << std::setw( 15 ) << std::setprecision( 6 ) <<
//    #define NL << '\n'
//    std::cout << "Original data:\n";
//    for ( size_t i = 0; i < xData.size(); i++ ) std::cout SP xData[i] SP yData[i] NL;
//    std::cout << "\nInterpolated data:\n";
//    for ( size_t i = 0; i < xVals.size(); i++ ) std::cout SP xVals[i] SP yVals[i] NL;

    Rocket sim("rocket.txt");

    return 0;
}