#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

#include "../RndGen/random.h"

using namespace std;

// Trial Wave Function
inline auto psi_T2 = [](double x, double mu, double sigma) {
    double esp1{(x - mu) / (sigma)};
    double esp2{(x + mu) / (sigma)};
    return pow(exp(-esp1 * esp1 / 2.) + exp(-esp2 * esp2 / 2.), 2);
};

// Kinetic divided by the trial wave function
inline auto Kin = [](double x, double mu, double sigma) {
    double xmu{x * mu};
    double sigma2{sigma * sigma};
    return -0.5 / (sigma2 * sigma2) *
           (x * x + mu * mu - sigma2 - 2 * xmu * tanh(xmu / sigma2));
};

// Potential Energy
inline auto Pot = [](double x) {
    double x2{x * x};
    return (x2 * x2 - 2.5 * x2);
};

// Hamiltonian divided by the trial wave function
inline auto Hamiltonian = [](double x, double mu, double sigma) {
    return Kin(x, mu, sigma) + Pot(x);
};

class VMC_simu
{
  private:
    Random rnd;                   // random number generator
    int _naccept{}, _nattempts{}; // for acceptance rate
    double _x, _delta, _delta_mu, _delta_sigma, _mu,
        _sigma; // starting point,step size and parameters of the trial wave
                // function
    double _mu_error{}, _sigma_error{}; // parameters of the trial wave function
    double _nblk{}, _nsteps{}, _steps{}, _block_ave{},
        _block_err{}; // blocks, total steps, steps in one block, block average
                      // and block error
    double _Lnew{}, _Lold{}, _Lnew_err{},
        _Lold_err{}; // cost function for simulated annealing
    double _T{};     // temperature for simulated annealing
    double _L_min{}, _L_min_err{}, _mu_min{},
        _sigma_min{}; // minimum value of the cost function and the
                      // corresponding parameters
  public:
    // Constructor: <starting point, step size, mu, sigma> and initialize the
    // random number generator
    VMC_simu();
    ~VMC_simu() = default;

    // input
    void Input();

    // Data blocking algorithm
    void Data_blocking(int blocks, int steps);

    // Metropolis algorithm: returns the new value of x
    double Metro(std::function<double(double, double, double)> f_psi2,
                 double x0, double delta);

    // Integrate a function using metropolis algorithm to sample the f_psi
    double Integrade(std::function<double(double, double, double)> f_psi2,
                     std::function<double(double, double, double)> fH,
                     int nsteps);

    // Simulated Annealing algorithm: return the value of the cost function
    // given a pdf and a cost function
    double Simu_annealing(double T);

    // Error function
    inline double Error(double ave, double ave2, int n) const
    {
        return n == 1 ? 0 : sqrt((ave2 - ave * ave) / (n - 1));
    }

    // Acceptance function
    inline double getacceptance() const
    {
        return static_cast<double>(_naccept) / _nattempts;
    }

    // Access operator var
    double &operator[](const std::string &var)
    {
        if (var == "nblk")
            return _nblk; // return the value of n block
        else if (var == "nsteps")
            return _nsteps; // return the value of n steps
        else if (var == "mu")
            return _mu; // return the value of mu
        else if (var == "sigma")
            return _sigma; // return the value sigma
        else if (var == "mu_error")
            return _mu_error; // return the value of mu error
        else if (var == "sigma_error")
            return _sigma_error; // return the value of sigma error
        else if (var == "block_ave")
            return _block_ave; // return the value of block ave
        else if (var == "block_err")
            return _block_err; // return the value of block err
        else if (var == "T")
            return _T; // return the value of Temperature
        else if (var == "delta_mu")
            return _delta_mu; // return the value of delta mu
        else if (var == "delta_sigma")
            return _delta_sigma; // return the value of delta sigma
        else if (var == "block_err")
            return _block_err; // return the value of block err
        else if (var == "L_min")
            return _L_min; // return the value of L min
        else if (var == "L_min_err")
            return _L_min_err; // return the value of L min err
        else if (var == "mu_min")
            return _mu_min; // return the value of mu min
        else if (var == "sigma_min")
            return _sigma_min; // return the value sigma m in
        // else if (var == "x")
        //    return _x; // return the value of x
        // else if (var == "Lnew")
        //    return _Lnew; // return the value of Lnew
        // else if (var == "Lold")
        //    return _Lold; // return the value of Lold
        // else if (var == "Lnew_err")
        //    return _Lnew_err; // return the value of Lnew_err
        // else if (var == "Lold_err")
        //    return _Lold_err; // return the value of Lold_err
        else
            throw std::invalid_argument(var + ": Invalid variable name");
    }
};