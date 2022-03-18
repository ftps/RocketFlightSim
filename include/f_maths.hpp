#ifndef F_MATH_HPP
#define F_MATH_HPP

#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <functional>
#include <iostream>

namespace fm {

    template<typename T>
    T sq(const T& v)
    {
        return v*v;
    }

    /*
    Our Vector data will simply be a relabelled std::array
    Special vectors include Vector3D and Vector2D, which represent
    vectors in R3 and R2 respectively
    */

    template<typename T, std::size_t N>
    using Vector = std::array<T,N>;
    using Vector3D = Vector<double,3>;
    using Vector2D = Vector<double,2>;

    // |v|²_2
    template<typename T, std::size_t N>
    T mod2V(const Vector<T,N>& v)
    {
        T res = 0;
        for(T vv : v) res += vv*vv;
        return res;
    }

    // |v|_2
    template<typename T, std::size_t N>
    double modV(const Vector<T,N>& v)
    {
        return sqrt(mod2V(v));
    }

    // |v|_p
    template<typename T, std::size_t N>
    double modVp(const Vector<T,N>& v, const double& p = 2)
    {
        double res = 0;
        for(T vv : v) res += std::pow(std::abs(vv), p);
        return std::pow(res, 1/p);
    }

    // dot v*v
    template<typename T, std::size_t N>
    inline double dot(const Vector<T,N>& l, const Vector<T,N>& r)
    {
        return l*r;
    }

    // cross vxv
    Vector3D cross(const Vector3D& l, const Vector3D& r)
    {
        Vector3D res;

        for(int i = 0; i < 3; ++i){
            res[i] = l[(i+1)%3]*r[(i+2)%3] - r[(i+1)%3]*l[(i+2)%3];
        }

        return res;
    }




    /*
    N1xN2 matrix class implementation
    Simple case of 
    */

    template<std::size_t N1, std::size_t N2, typename T>
    class Matrix{
        private:
            Vector<Vector<T, N2>, N1> v;
        public:
            Matrix() : v({0}) { }
            Matrix(const Vector<Vector<T,N2>,N1>& v) : v(v) { }
            Vector<T, N2> operator[](const int& n) const { return v[n]; }
            Vector<T, N2>& operator[](const int& n) { return v[n]; }

            Matrix<N1,N2,T> operator+(const Matrix<N1,N2,T>& r)
            {
                Matrix<N1,N2,T> res;
                for(std::size_t i = 0; i < N1; ++i){
                    for(std::size_t j = 0; j < N2; ++j){
                        res[i][j] = v[i][j]+r[i][j];
                    }
                }
                return res;
            }

            Matrix<N1,N2,T> operator-(const Matrix<N1,N2,T>& r)
            {
                Matrix<N1,N2,T> res;
                for(std::size_t i = 0; i < N1; ++i){
                    for(std::size_t j = 0; j < N2; ++j){
                        res[i][j] = v[i][j]-r[i][j];
                    }
                }
                return res;
            }
    };

    using Matrix22 = Matrix<2,2,double>;
    using Matrix33 = Matrix<3,3,double>;

    template<std::size_t N1, std::size_t N2, typename T>
    double modMpq(const Matrix<N1,N2,T>& M, const double& p, const double& q)
    {
        double res = 0, aux;

        for(std::size_t i = 0; i < N1; ++i){
            aux = 0;
            for(std::size_t j = 0; j < N2; ++j){
                aux += std::pow(std::abs(M[i][j]), p);
            }
            res += std::pow(aux, q/p);
        }
        return std::pow(res, 1/q);
    }
};

/*
Vector print operator
*/
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const fm::Vector<T,N>& v)
{
    os << "[ ";
    for(T vv : v){
        os << vv << " ";
    }
    os << "]";
    return os;
}

/*

Vector mathematical operators:
vector * vector (dot product)
vector + vector, +=
vector - vector, -= and neg
scalar*vector, vector*scalar (same)
vector/scalar
2-norm squared, 2-norm, p-norm
cross and dot functions

*/

// dot product
template<typename T, std::size_t N>
T operator*(const fm::Vector<T,N>& l, const fm::Vector<T,N>& r)
{
    T res = 0;

    for(std::size_t i = 0; i < N; ++i){
        res += l[i]*r[i];
    }

    return res;
}

// + operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator+(const fm::Vector<T,N>& l, const fm::Vector<T,N>& r)
{
    fm::Vector<T,N> res;

    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i] + r[i];
    }

    return res;
}

// += operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator+=(fm::Vector<T,N>& l, const fm::Vector<T,N>& r)
{
    for(std::size_t i = 0; i < N; ++i){
        l[i] += r[i];
    }

    return l;
}

// - operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator-(const fm::Vector<T,N>& l, const fm::Vector<T,N>& r)
{
    fm::Vector<T,N> res;

    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i] - r[i];
    }

    return res;
}

// negation operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator-(const fm::Vector<T,N>& v)
{
    fm::Vector<T,N> res;

    for(std::size_t i = 0; i < N; ++i){
        res[i] = -v[i];
    }

    return res;
}

// -= operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator-=(fm::Vector<T,N>& l, const fm::Vector<T,N>& r)
{
    for(std::size_t i = 0; i < N; ++i){
        l[i] -= r[i];
    }

    return l;
}

// v*scalar operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator*(const fm::Vector<T,N>& l, const T& r)
{
    fm::Vector<T,N> res;

    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i]*r;
    }

    return res;
}

// scalar*v operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator*(const T& l, const fm::Vector<T,N>& r)
{
    return r*l;
}

// v/scalar operator
template<typename T, std::size_t N>
fm::Vector<T,N> operator/(const fm::Vector<T,N>& l, const T& r)
{
    fm::Vector<T,N> res;

    for(std::size_t i = 0; i < N; ++i){
        res[i] = l[i]/r;
    }

    return res;
}

template<std::size_t N1, std::size_t N2, typename T>
std::ostream& operator<<(std::ostream& os, const fm::Matrix<N1,N2,T>& m)
{
    os << "[";
    for(std::size_t i = 0; i < N1; ++i){
        for(std::size_t j = 0; j < N2; ++j){
            os << " " << m[i][j];
        }
        os << ((i == N1-1) ? " ]" : "\n ");
    }
    return os;
}

template<std::size_t N1, std::size_t N2, typename T>
fm::Matrix<N1,N2,T> operator+=(fm::Matrix<N1,N2,T>& l, const fm::Matrix<N1,N2,T>& r)
{
    for(std::size_t i = 0; i < N1; ++i){
        l[i] += r[i];
    }

    return l;
}

template<std::size_t N1, std::size_t N2, typename T>
fm::Matrix<N1,N2,T> operator-=(fm::Matrix<N1,N2,T>& l, const fm::Matrix<N1,N2,T>& r)
{
    for(std::size_t i = 0; i < N1; ++i){
        l[i] -= r[i];
    }

    return l;
}

template<std::size_t N1, std::size_t N2, std::size_t N3, typename T>
fm::Matrix<N1,N2,T> operator*(const fm::Matrix<N1,N3,T>& l, const fm::Matrix<N3,N2,T>& r)
{
    fm::Matrix<N1,N2,T> res;

    for(std::size_t i = 0; i < N1; ++i){
        for(std::size_t j = 0; j < N2; ++j){
            for(std::size_t k = 0; k < N3; ++k){
                res[i][j] += l[i][k]*r[k][j];
            }
        }
    }

    return res;
}


template<std::size_t N1, std::size_t N2, typename T>
fm::Vector<T,N1> operator*(const fm::Matrix<N1,N2,T>& l, const fm::Vector<T,N2>& r)
{
    fm::Vector<T,N2> res;

    for(std::size_t i = 0; i < N1; ++i){
        res[i] = l[i]*r;
    }

    return res;
}

template<std::size_t N1, std::size_t N2, typename T>
fm::Vector<T,N2> operator*(const fm::Vector<T,N1>& l, const fm::Matrix<N1,N2,T>& r)
{
    fm::Vector<T,N2> res;

    for(std::size_t i = 0; i < N2; ++i){
        for(std::size_t j = 0; j < N1; ++j){
            res[i] = l[j]*r[j][i];
        }
    }

    return res;
}

namespace fm {

    template<typename T>
    using pair_tV = std::vector<std::pair<double, T>>; // vector of pairs of a double and another variable (t, V)
    using XY = std::pair<double, double>;   // point x,y
    using pairXY = std::vector<XY>; // vectors of points x,y (easy for gnuplot)

    template<typename T>
    using function_tV = std::function<T(const double&, const T&)>;

    template<typename T>
    using function_V = std::function<T(const T&)>;

    template<typename T>
    using endCond = std::function<const bool(const double&, const T&)>;


    using funcInt_data = std::function<double(const pairXY&)>;
    using func_RR = std::function<double(const double&)>;
    using funcInt_func = std::function<double(const func_RR&, const XY&, const int&)>;

    /*
    Integrators for ODEs
    */
    
    namespace ode {

        #define STEP_PARAMS const function_tV<T>& f, const double& t, const T& y, const double& dt

        template<typename T>
        using ODE_step = std::function<T(const function_tV<T>&, const double&, const T&, const double&)>;

        template<typename T>
        T EulerExp(STEP_PARAMS)
        {
            return y + f(t, y)*dt;
        }

        template<typename T>
        T RK2(STEP_PARAMS)
        {
            T k1, k2;

            k1 = f(t,y)*dt;
            k2 = f(t+dt, y + k1*dt);

            return y + 0.5*(k1 + k2); 
        }

        template<typename T>
        T RK4(STEP_PARAMS)
        {
            T k1, k2, k3, k4;
            double dt2 = dt/2;

            k1 = f(t,y)*dt;
            k2 = f(t+dt2, y + k1*dt2);
            k3 = f(t+dt2, y + k2*dt2);
            k4 = f(t+dt, y + k3*dt);

            return y + (k1+2*k2+2*k3+k4)*(dt/6);
        }

        /*
        ODE solver options

        dt - time-step
        save_step - saves data in step intervals of this size
        endCond - end conditions function, return true to continue and false to stop;
        step - type of solver step (Euler Explicite, Runge-Kutta 4, etc)
        */
        template<typename T>
        struct ODEsolve_opt {
            double dt = 1e-3;
            uint64_t save_step = 1;
            endCond<T> end = [](const double& t, const T&){ return t <= 10; };
            ODE_step<T> step = EulerExp<T>;
        };

        // main ODE solver
        template<typename T>
        pair_tV<T> ODEsolve(const function_tV<T>& f, const double& t0, const T& y0, const ODEsolve_opt<T>& opts = ODEsolve_opt<T>())
        {
            pair_tV<T> res;            
            double t = t0;
            T y = y0;
            uint64_t save_count = 0;

            std::cout << "Starting integration . . ." << std::endl;
            do{
                if(!save_count){
                    res.emplace_back(t, y);
                }

                y = opts.step(f, t, y, opts.dt);

                save_count = (save_count+1)%opts.save_step;
                t += opts.dt;
            }while(opts.end(res.back().first, res.back().second) || save_count);
           
            std::cout << "Done." << std::endl;

            return res;
        }

        

        /*
        Sympletic integrators
        Two equal variables, one representing speed and another position
        Function for acceleration as a function of time and position
        */

        template<typename T>
        void Leapfrog(T& x, T& v, const function_V<T>& a, const double& dt)
        {
            T ax = a(x);
            x += v*dt + ax*(0.5*sq(dt));
            v += (ax + a(x))*(0.5*dt);
        }

        //template<typename T>
        //using function_V = std::function<T(const T&)>;
        // std::vector<double>
        template<typename T, typename U>
        void Leapfrog(T& x, T& v, const U& a, const double& dt, T& a2)
        {
            T ax = a2;
            x += v*dt + ax*(0.5*sq(dt));
            a2 = a(x);
            v += (ax + a2)*(0.5*dt);
        }



    };

    /* Regular integrators for data and functions */
    namespace intgr {

        /*
        Forward Riemann Sum integration
        */

        double RS_ForD(const pairXY& data)
        {
            double res = 0;

            for(std::size_t i = 0; i < data.size()-1; ++i){
                res += data[i+1].second*(data[i+1].first - data[i].first);
            }

            return res;
        }

        double RS_ForF(const func_RR& f, const XY& ep, const int& N)
        {
            double dt = (ep.second - ep.first)/(N-1);
            double res = 0;

            for(double t = ep.first; (ep.second > ep.first) ? t < ep.second : t > ep.second; t += dt){
                res += f(t+dt);
            }

            return res*dt;
        }


        /*
        Backwards Riemann Sum integration
        */

        double RS_BackD(const pairXY& data)
        {
            double res = 0;

            for(std::size_t i = 0; i < data.size()-1; ++i){
                res += data[i].second*(data[i+1].first - data[i].first);
            }

            return res;
        }

        double RS_BackF(const func_RR& f, const XY& ep, const int& N)
        {
            double dt = (ep.second - ep.first)/(N-1);
            double res = 0;

            for(double t = ep.first; (ep.second > ep.first) ? t < ep.second : t > ep.second; t += dt){
                res += f(t);
            }

            return res*dt;
        }


        /*
        Trapezoid Rule integration
        */

        double TrapD(const pairXY& data)
        {
            double res = 0;

            for(std::size_t i = 0; i < data.size()-1; ++i){
                res += 0.5*(data[i+1].second + data[i].second)*(data[i+1].first - data[i].first);
            }

            return res;
        }

        double TrapF(const func_RR& f, const XY& ep, const int& N)
        {
            double dt = (ep.second - ep.first)/(N-1);
            double res = 0.5*f(ep.first);

            for(double t = ep.first + dt; (ep.second > ep.first) ? t <= (ep.second-dt) : t >= (ep.second-dt); t += dt){
                res += f(t);
            }

            return (res + 0.5*f(ep.second))*dt;
        }


        /*
        Simpson's rule integration
        */

        double SimpsonD(const pairXY& data)
        {
            double res = 0;
            double h2, h21;

            for(std::size_t i = 0; i < ((data.size()-1)/2 - 1); ++i){
                h2 = data[2*i+1].first - data[2*i].first;
                h21 = data[2*(i+1)].first - data[2*(i+1)-1].first;
                res += ((h21 + h2)/6)*((2 - h21/h2)*data[2*i].second + (sq(h21 + h2)/(h21*h2))*data[2*i+1].second + (2 - h2/h21)*data[2*(i+1)].second);
            }

            if(!(data.size()%2)){
                double a, b, c;
                std::size_t i = data.size()-2;

                h2 = data[i].first - data[i-1].first;
                h21 = data[i+1].first - data[i].first;
            
                a = (2*sq(h21) - 3*h21*h2)/(6*(h21+h2));
                b = (sq(h21) + 3*h21*h2)/(6*h2);
                c = h21*sq(h21)/(6*h2*(h21 +h2));

                res += a*data[i+1].second + b*data[i].second + c*data[i-1].second;
            }

            return res;
        }

        double SimpsonF(const func_RR& f, const XY& ep, const int& N)
        {
            double dt = (ep.second - ep.first)/(N-1);
            double res = 0.5*f(ep.first);

            for(double t = ep.first + dt; (ep.second > ep.first) ? t <= (ep.second-dt) : t >= (ep.second-dt); t += dt){
                res += f(t);
            }

            return (res + 0.5*f(ep.second))*dt;
        }
        

        /*
        Main integrator functions
        */

        double integrate(const pairXY& data, const funcInt_data& f = TrapD)
        {
            return f(data);
        }

        double integrate(const func_RR& f, const XY& ep, const int& N, const funcInt_func& ff = TrapF)
        {
            return ff(f, ep, N);
        }

    };
};

#endif