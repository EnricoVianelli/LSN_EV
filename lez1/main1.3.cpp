
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

  //------------test chi--------------------------------------
  double chi = 0;// statistical chi-sqrd 
  double chi_j;// chi-sqrd on single interval
  double x;// rnd number
  int M = 100;   //# of intervals 
  int N = 10000; //# of throws
  int exp = N / M;
  double L = 1. / double(M);// with of sub-intervals
  int counter;// acc
  ofstream WriteData;
  WriteData.open("OUTPUT/Test-Chi.out");

  for (int j = 0; j < M; j++) {
    counter = 0;
    chi_j = 0;
    for (int i = 0; i < N; i++) {// throw 10000 numbers btw 0 and 1
      
      x = rnd.Rannyu();
        if (x >= L * j && x < L * (j + 1)) {
        counter++;
        }
    }
    chi_j += double(counter-exp)*double(counter-exp) / double(exp);
    chi += chi_j;

    // ....saving data....
    
    WriteData << j+1 << " " << chi_j << endl;
  }

  cout << "Chi-square:" << chi << endl;

  WriteData.close();
  rnd.SaveSeed();
  return 0;
}
