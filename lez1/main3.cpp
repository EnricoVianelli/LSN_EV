#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "random.h"  // Your Random class header

using namespace std;

int main(int argc, char *argv[]) {

    Random rnd;
    int seed[4];
    int p1, p2;

    // 1) Read primes from file
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();

    // 2) Initialize random seed from file
    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    // --- Buffon's experiment setup ---
    double d = 1.0;    // distance between lines
    double L = 0.4;    // needle length
    int N = 100;       // number of blocks
    int M = 1000000;  // total throws
    int A = M / N;     // throws per block

    int M_hit = 0, M_miss = 0;
    int M_throws = 0;

    double pi_ave = 0.0, pi_ave2 = 0.0;
    double stat_error = 0.0;

    // Prepare output
    ofstream Wdata("OUTPUT/pi_"+to_string(A)+"perblock.out");
    Wdata << setw(12) << "# BLOCK:" 
          << setw(16) << "PI:" 
          << setw(16) << "ERROR:" << endl; 

    // --- Main loop over N blocks ---
    for(int j = 0; j < N; j++) {

        // Perform A throws in block j
        for(int i = 0; i < A; i++) {

            // 1) Sample the center of the needle 
            double x = rnd.Rannyu(0, d*0.5);

            // 2) Sample angle theta in [0, π), then compute half-length projection
            double theta = rnd.Rannyu(0, M_PI);
            double r  = (L * 0.5) * sin(theta);

            // 3) Check if needle crosses line
            if (x <= r) M_hit++;
            else        M_miss++;

            M_throws++;
        }

        // Probability of hitting
        double P = double(M_hit) / double((j+1) * A);

        // Buffon's formula
        double pi_temp = (2.0 * L) / (d * P);

        // Accumulate for averages
        pi_ave  += pi_temp;
        pi_ave2 += pi_temp * pi_temp;

        // Compute statistical error
        double ave  = pi_ave  / double(j+1);
        double ave2 = pi_ave2 / double(j+1);

        if (j == 0) {
            // For the 1st block, you can't get a meaningful error
            stat_error = 0.0;
        } else {
            stat_error = sqrt((ave2 - ave*ave) / double(j));
        }

        // Write out results
        Wdata << setw(12) << j+1
              << setw(16) << ave
              << setw(16) << stat_error << endl;
    }

    Wdata.close();

    // Quick check
    cout << "hit:    " << M_hit << endl;
    cout << "miss:   " << M_miss << endl;
    cout << "throws: " << M_throws << endl;

    // Save final seed
    rnd.SaveSeed();

    return 0;
}
