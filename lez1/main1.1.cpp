
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

  //------------INTEGRAL--------------------------------------------
  double x;
  int M = 100000;
  int N = 100;
  int L = M / N;
  double sum;
  double ave;
  double sum_prog = 0;
  double sum_progquad = 0;
  double avequad;
  double staterror;

  ofstream WriteData;
  WriteData.open("OUTPUT/Integrali-1.1.out");

  for (int j = 0; j < N ; j++) {

    sum = 0;
    ave = 0;
    avequad = 0;
    staterror = 0;

    for (int i = 0; i < L; i++) {

      x = rnd.Rannyu();
      sum += x;
    }
    ave = sum / double(L); // Compute average for the current block
    avequad = ave * ave;   // Compute average squares
    sum_prog += ave; // sum of averages
    sum_progquad += avequad; // sum of averages squares

    // Compute statistical error
    if(j==0){
      staterror = 0.0;
    }else{
      staterror = sqrt((sum_progquad / double(j + 1) -
                        (sum_prog / double(j + 1)) * (sum_prog / double(j + 1))) /double(j));
    }

    WriteData << j+1 << " " << sum_prog / (j + 1) << " " << staterror << endl;
  }
  WriteData.close();

  //saving seed 
  rnd.SaveSeed();
  return 0;
}
