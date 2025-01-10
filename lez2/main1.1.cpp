
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

double function(double x) { return (M_PI / 2) * cos(M_PI * x / 2); }

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
  int M = 10000;
  int N = 100;
  int L = M / N;
  double sum, ave , avequad; // tmp accumolator
  double sum_prog = 0;
  double sum_progquad = 0;
  double staterror , mean , mean2;// mean over blocks

  ofstream WriteData;
  WriteData.open("Integrali-1.out");

  for (int j = 0; j < N; j++) {

    sum = 0;
    ave = 0;
    avequad = 0;
    staterror = 0;

    for (int i = 0; i < L; i++) {

      x = rnd.Rannyu();
      sum += function(x);
    }
    ave = sum / double(L); // Compute average for the current block
    avequad = ave * ave;   // Compute average of squares
    sum_prog += ave;
    sum_progquad += avequad;
    // statistical errorr
    mean = sum_prog/double(j+1);
    mean2 = sum_progquad/double(j+1);
    if(j==0){
      staterror = 0.0;
      //sqrt((sum_progquad / double(j + 1) -
                        //(sum_prog / double(j + 1)) * (sum_prog / double(j + 1))) /1.);
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
