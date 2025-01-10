
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
  //-------------------Compute market price -----------
  //problem param
  double T = 1.;
  double K = 100.;
  double r = 0.1;
  double sigma = 0.25;
  double t0, S, S0, t_current ;
  //data blocking variables
  double c_ave, p_ave, c_ave2, p_ave2, c_staterror, p_staterror;
  double c_sprog = 0;
  double p_sprog = 0;
  double c_sprog2 = 0;
  double p_sprog2 = 0;
  int N =100;// # of blocks
  int A = 2000 ;// throws per block

  ofstream callout, putout;//printing results
  callout.open("Discr_call.out");
  putout.open("Discr_put.out");

  double t_i = T/100.;
//---------------loop over blocks---------------------------
  for (int j = 0; j < N; j++) {
    c_ave = 0.;
    p_ave = 0.;
    p_ave2 = 0.;
    c_ave2 = 0.;
    c_staterror = 0.;
    p_staterror = 0.;
    // discrete sampling 
    for (int i = 0; i < A; i++) {
        S0 = 100.;
        t0 = 0.;
        t_current = t_i;
        // discrete sampling
        for(int k = 0; k<100; k++){//capire perchè non va ?
            S = S0*exp((r-0.5*sigma*sigma)*(t_current-t0)+sigma*sqrt(t_current-t0)*rnd.Gauss(0.,1.));
            S0 = S;
            t0 += t_i;
            t_current += t_i; 
            //if (k==99) cout <<"t_c:" << t_current << "\trestart" <<  endl;
        }
        c_ave += exp(-r*T)*findMax(0,S-K);
        p_ave += exp(-r*T)*findMax(0,K-S);
    }
    // Compute average for the current block
    c_ave = c_ave / double(A); 
    p_ave = p_ave / double(A); 
    // Compute average sqr for the current block
    c_ave2 = c_ave * c_ave;   
    p_ave2 = p_ave * p_ave;
    // Blocks avg
    c_sprog += c_ave;
    p_sprog += p_ave;
    c_sprog2 += c_ave2;
    p_sprog2 += p_ave2;
    // statistical error
    if(j==0){
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
    putout << scientific << p_sprog/double(j+1) << " " << p_staterror << endl;
  }
  cout << "Blocks:" << N << "\tTotal throws:" << A*N << "\tBlocks size:" << A << endl;
  callout.close();
  putout.close();
  rnd.SaveSeed();

  return 0;
}