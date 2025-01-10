#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  Random rnd;
  int seed[4];
  int p1, p2;

  ifstream Primes("Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

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
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;

  //------------N=1------------------------------------
  double x_e, x_d, x_l;
  int M = 10000;
  int N[4] = {1, 2, 10, 100};
  ofstream Dout, Eout, Lout;
  for(int k = 0; k < 4; k++){

    Dout.open("OUTPUT/Dice_distribution" + std::to_string(N[k]) + ".out");
    Eout.open("OUTPUT/Exp_distribution" + std::to_string(N[k]) + ".out");
    Lout.open("OUTPUT/Lorenyz_distribution" + std::to_string(N[k]) + ".out");
    for (int j = 0; j < M; j++){
    // test expon and lorentz
      x_d = 0;
      x_e = 0;
      x_l = 0;

      for(int i = 0; i < N[k]; i++){
        x_d += rnd.Rannyu();
        x_e += rnd.Expo(1.);
        x_l += rnd.Lorentz(0., 1.);
      }

      Dout << x_d/N[k] << endl;
      Eout << x_e/N[k] << endl;
      Lout << x_l/N[k] << endl;
    }

    Dout.close();
    Eout.close();
    Dout.close();
  }

  rnd.SaveSeed();// saving seed
  return 0;
}
