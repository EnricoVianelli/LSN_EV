#include "VMC_simu.h"

VMC_simu::VMC_simu()
{
    Input();
    rnd.initialize();
    // initialize the cost function for simulated annealing
    Data_blocking(_nblk, _steps);
    _Lnew = _block_ave;
    _Lnew_err = _block_err;
    // initialize the minimum value of the cost function
    _L_min = _Lnew;
    _mu_min = _mu;
    _sigma_min = _sigma;
}

void VMC_simu::Input()
{
    std::ifstream Rinput;

    // Read input informations
    Rinput.open("input.dat");

    if (!Rinput.is_open())
    {
        std::cerr << "Error: unable to open input.dat\n";
        std::exit(1);
    }

    std::string line;
    std::string property;
    while (getline(Rinput, line))
    {
        stringstream ss(line); // Usa stringstream per analizzare la riga
        ss >> property;        // Leggi la prima parola della riga
        if (property == "N_BLK:")
            ss >> _nblk;
        else if (property == "N_STEP_TOT:")
            ss >> _nsteps;
        else if (property == "X_0:")
            ss >> _x;
        else if (property == "DELTA_X:")
            ss >> _delta;
        else if (property == "MU:")
            ss >> _mu;
        else if (property == "SIGMA:")
            ss >> _sigma;
        else if (property == "T_0:")
            ss >> _T;
        else if (property == "DELTA_MU:")
            ss >> _delta_mu;
        else if (property == "DELTA_SIGMA:")
            ss >> _delta_sigma;
        else if (property == "MU_ERROR:")
            ss >> _mu_error;
        else if (property == "SIGMA_ERROR:")
            ss >> _sigma_error;
        else
            cerr << line << " unknown input parameter" << endl;
    }
    _steps = _nsteps / _nblk;

    std::cout << "Total number of steps = " << _nsteps << '\n';
    std::cout << "Number of blocks = " << _nblk << '\n';
    std::cout << "Number of steps in one block = " << _steps << '\n' << '\n';
    std::cout << "Step lenght = " << _delta << '\n';
    std::cout << "Initial position = " << _x << '\n';
    std::cout << "Initial mu = " << _mu << '\n';
    std::cout << "Initial sigma = " << _sigma << '\n';
    std::cout << "Initial temperature = " << _T << '\n';
    std::cout << "Step size for mu = " << _delta_mu << '\n';
    std::cout << "Step size for sigma = " << _delta_sigma << '\n';
    std::cout << "Desired error for mu = " << _mu_error << '\n';
    std::cout << "Desired error for sigma = " << _sigma_error << '\n';

    Rinput.close();
}
double VMC_simu::Simu_annealing(double T)
{

    _Lold = _Lnew;
    _Lold_err = _Lnew_err;
    _T = T;

    double mu_old{_mu};
    double sigma_old{_sigma};

    _mu += rnd.Rannyu(-_delta_mu, _delta_mu) * T;

    do
    {
        _sigma += rnd.Rannyu(-_delta_sigma, _delta_sigma) * T;

    } while (fabs(_sigma) < 0.002); // sigma too small the wave function
                                    // diverges

    Data_blocking(_nblk, _steps);
    _Lnew = _block_ave;
    _Lnew_err = _block_err;

    if (_Lnew < _L_min)
    {
        _L_min = _Lnew;
        _L_min_err = _Lnew_err;
        _mu_min = _mu;
        _sigma_min = _sigma;
    }

    // Metropolis algorithm for SA
    double acceptance = exp(-(_Lnew - _Lold) / T);

    if (acceptance >= 1) { return _Lnew; }
    else if (rnd.Rannyu() < acceptance)
    {
        return _Lnew;
    }

    _mu = mu_old;
    _sigma = sigma_old;
    return _Lold;
}

/**
  metropolis for sampling the trial wave function
  **/
double VMC_simu::Metro(std::function<double(double, double, double)> f_psi2,
                       double x0, double delta)
{

    _nattempts++;
    double x = x0 + rnd.Rannyu(-delta, delta);
    double p = f_psi2(x, _mu, _sigma) / f_psi2(x0, _mu, _sigma);
    if (p >= 1)
    {
        _naccept++;
        return x;
    }
    else
    {
        if (rnd.Rannyu() < p)
        {
            _naccept++;
            return x;
        }
    }
    return x0;
}

void VMC_simu::Data_blocking(int nblocks, int nsteps)
{
    double ave{}, sum{}, sum2{};
    for (int i{}; i < nblocks; i++)
    {
        ave = Integrade(psi_T2, Hamiltonian, nsteps);
        sum += ave;
        sum2 += ave * ave;
    }
    _block_ave = sum / nblocks;
    _block_err = Error(sum / nblocks, sum2 / nblocks, nblocks);
}

/**
 Integrate over nsteps of one blk sampling square mod of WF
 **/
double VMC_simu::Integrade(std::function<double(double, double, double)> f_psi2,
                           std::function<double(double, double, double)> fH,
                           int nsteps)
{
    double integral{};
    for (int i{}; i < nsteps; i++)
    {

        _x = this->Metro(f_psi2, _x, _delta);
        integral += fH(_x, _mu, _sigma);
    }

    return integral / nsteps;
}
