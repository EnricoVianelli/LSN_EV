
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include "VMC_simu.h"

using namespace std;

int main(int argc, char *argv[])
{

    // SA
    array<string, 2> cooling_laws = {"exponential", "linear"};

    if (argc != 2 || atoi(argv[1]) < 0 || atoi(argv[1]) > 1)
    {

        if (argc != 2)
            cerr << "Usage: " << argv[0] << "  <code cooling law>" << endl;
        else
            cerr << " Invalid code cooling law !" << endl;

        cerr << " Available code:" << endl;

        for (int i = 0; i < (int)std::size(cooling_laws); i++)
            cerr << "\t" << i << " = " << cooling_laws[i] << endl;

        cerr << endl;
        return 1;
    }

    int cooling = atoi(argv[1]);
    VMC_simu simu{};

    cout << "Cooling law = " << cooling_laws[cooling] << "\n";

    fstream Tout{"OUTPUT/" + cooling_laws[cooling] + "/temperature.dat",
                 ios::out};
    fstream Hout{"OUTPUT/" + cooling_laws[cooling] + "/energy.dat", ios::out};
    fstream paramsout{"OUTPUT/" + cooling_laws[cooling] + "/parameters.dat",
                      ios::out};

    int k{}; // SA iteration
    double T{simu["T"]};
    const double delta_mu{simu["delta_mu"]};
    const double delta_sigma{simu["delta_sigma"]};
    const double err_mu{simu["mu_error"]};
    const double err_sigma{simu["sigma_error"]};

    // When the step size in mu or sigma is smaller than the desired error then
    // stop because T is too low
    Tout << setw(12) << "Iteration:" << setw(12) << "T:\n";
    Hout << setw(12) << "Iteration:" << setw(16) << "H(SA):" << setw(16)
         << "error:\n";
    paramsout << setw(16) << "mu:" << setw(16) << "sigma" << setw(16)
              << "energy:\n";
    while (2 * delta_mu * T > err_mu or 2 * delta_sigma * T > err_sigma)
    {
        switch (cooling)
        {
        case 0:
            T = T * std::pow(0.9999, k);
            break;
        case 1:
            T = T - 1e-5 * k;
            break;
        default:
            cerr << "Invalid cooling law\n";
            return 1;
        }

        Tout << setw(12) << k << setw(12) << T << "\n";

        simu.Simu_annealing(T);

        Hout << setw(12) << k << setw(16) << simu["block_ave"] << setw(16)
             << simu["block_err"] << "\n";
        paramsout << setw(16) << simu["mu"] << setw(16) << simu["sigma"]
                  << setw(16) << simu["block_ave"] << "\n";

        k++;
    }

    cout << "\nCompleted Simulated Annealing\t"
         << "\n";
    cout << "Minimum H = " << simu["L_min"] << "\nmu = " << fabs(simu["mu_min"])
         << "\nsigma = " << fabs(simu["sigma_min"]) << "\n";

    double H{};
    double sum_H{}, sum2_H{};

    simu["mu"] = simu["mu_min"];
    simu["sigma"] = simu["sigma_min"];
    // cout << simu.getmu() << "\n";
    ofstream wave{"OUTPUT/" + cooling_laws[cooling] + "/wavef_samples.dat",
                  ios::out};
    ofstream Hmin{"OUTPUT/" + cooling_laws[cooling] + "/H_min.dat", ios::out};
    cout << "\nSample the wave function with best parameters\n";
    wave << "X_samples of WF:\n";
    Hmin << setw(12) << "Block:" << setw(16) << "H_min:" << setw(16)
         << "error:\n";
    // Sampling of the modulus square of the wave function
    double x{};
    for (int i{1}; i <= simu["nsteps"] * 2; i++)
    {
        x = simu.Metro(psi_T2, x, 3);
        wave << x << "\n";
    }
    cout << "Starting to data block the value of the Energy\n";
    // Data blocking to get the value of the Energy
    for (int i{1}; i <= simu["nblk"]; i++)
    {
        H = simu.Integrade(psi_T2, Hamiltonian, 2000);
        sum_H += H;
        sum2_H += H * H;
        Hmin << setw(12) << i << setw(16) << sum_H / i << setw(16)
             << simu.Error(sum_H / i, sum2_H / i, i) << "\n";
    }

    wave.close();
    Hmin.close();
    sum_H /= simu["nblk"];
    sum2_H /= simu["nblk"];

    cout << '\n'
         << "\tH:\t" << sum_H << "\terror:\t"
         << simu.Error(sum_H, sum2_H, simu["nblk"]) << "\n";
    // cout << '\n' << "\tH:\t" << simu["block_ave"] << " +/- " <<
    // simu["block_err"] << "\n";

    return 0;
}
