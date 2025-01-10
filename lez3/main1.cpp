
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

double findMax(double a, double b) {
    return (a > b) ? a : b;
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

  //-------------------compute market price at time T-----------
  //problem param
  double T = 1.;    // time
  double S0 = 100.;  // starting asset preice    
  double K = 100.;   // strike price
  double r = 0.1;   // free risck interest rate
  double sigma = 0.25;  // volatility
  double t0 = 0.0;   //  starting time 
  double S ;        // asset prize variable
  //data blocking variables
  double c_ave, p_ave, c_ave2, p_ave2, c_staterror, p_staterror;
  double c_sprog = 0.;
  double p_sprog = 0.;
  double c_sprog2 = 0.;
  double p_sprog2 = 0.;
  int N =100;// # of blocks
  int A = 2000 ;// throws per block

  ofstream callout, putout;//printing results
  callout.open("Direct_call.out");
  putout.open("Direct_put.out");
//---------------main loop over blocks-----------------
  for (int j = 0; j < N; j++) {
    c_ave = 0.;
    p_ave = 0.;
    p_ave2 = 0.;
    c_ave2 = 0.;
    c_staterror = 0.;
    p_staterror = 0.;

    for (int i = 0; i < A; i++) {
      S = S0*exp((r-0.5*sigma*sigma)*(T-t0)+sigma*sqrt(T-t0)*rnd.Gauss(0,1));
      c_ave += exp(-r*T)*findMax(0,S-K);
      p_ave += exp(-r*T)*findMax(0,K-S);
    }
    // Compute average for the current block
    c_ave = c_ave / double(A); 
    p_ave = p_ave / double(A);
    // Compute squared average for the  current block
    c_ave2 = c_ave * c_ave;   
    p_ave2 = p_ave * p_ave;
    // avg over blocks 
    c_sprog += c_ave;
    p_sprog += p_ave;
    c_sprog2 += c_ave2;
    p_sprog2 += p_ave2;
    // statistical error
    if(j==0){// j<=1
      c_staterror = 0.0;
      p_staterror = 0.0;
    }else{
      c_staterror = sqrt((c_sprog2/ double(j + 1) -
                        (c_sprog / double(j + 1)) * (c_sprog / double(j + 1))) /double(j));
      p_staterror = sqrt((p_sprog2/ double(j + 1) -
                        (p_sprog / double(j + 1)) * (p_sprog / double(j + 1))) /double(j));
    }
    // ....saving data....
    callout << scientific << c_sprog/double(j+1) << " " << c_staterror << endl;
    putout  << scientific << p_sprog/double(j+1) << " " << p_staterror << endl;
  }
 cout << "Blocks:" << N << "\tTotal throws:" << A*N << "\tBlocks size:" << A << endl;
  callout.close();
  putout.close();
  rnd.SaveSeed();

  return 0;
}