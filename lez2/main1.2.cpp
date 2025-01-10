
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

double function(double x) { return M_PI *0.5 * cos(M_PI * x *0.5); }

double prob_function1(double x) {   return 2.*(1.-x); }

double prob_function2(double x) { // something like the first order expansion 
    double N = M_PI*(2-M_PI/48.);// normalizazion
    return (M_PI/N)*2.*(1.-(M_PI*x*x)/16.);
}

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
  double x, y;
  double weight;
  //double fmax = 2. ;
  int M = 10000;
  int N = 100;
  int L = M / N;
  double sum, ave, avequad;
  double sum_prog = 0;
  double sum_progquad = 0;
  double staterror, mean, mean2;// mean over the blocks

  ofstream WriteData;
  WriteData.open("Integrali-2.out");

  for (int j = 0; j < N; j++) {

    sum = 0;
    ave = 0;
    avequad = 0;
    staterror = 0;

    for (int i = 0; i < L; i++) {
        // Sampling x by Accept Reject
        /*do {
            x = rnd.Rannyu();
            y = rnd.Rannyu();
        }while(y < prob_function1(x)/fmax);*/
        // inverse function
        y = rnd.Rannyu();
        x = 1-sqrt(1-y);
        // Calculate the correction btw uniform distr and my distribution
        weight = 1./prob_function1(x);// slide better
        sum += weight*function(x);
    }
    ave = sum / double(L); // Compute average for the current block
    avequad = ave * ave;   // Compute average of squares
    sum_prog += ave;
    sum_progquad += avequad;

    // computing stat error
    mean = sum_prog/double(j+1);
    mean2 = sum_progquad/double(j+1);

    if(j==0){
      staterror = 0.0; 
    }else{
      staterror = sqrt((mean2 -mean*mean) /double(j));
    }
    // ....saving data....
    WriteData << j+1  << " " << mean << " " << staterror
              << endl;
  }
  WriteData.close();
  //saving seed
  rnd.SaveSeed();
  return 0;
}