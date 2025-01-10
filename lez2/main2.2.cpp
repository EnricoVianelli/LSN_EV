#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>

using namespace std;

double get_squared_modulus(vector <double> v) { 
    if(v.size()==3) {// check if v is a 3D coord vector
        double x = v[0];
        double y = v[1];
        double z = v[2];
        return (x*x + y*y + z*z);
    }else
    cerr << "Error: your vector is not a 3D vector" << endl;
    return 0;
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
  //------------Random Walk in continuom---------------------------
  double a = 1.; // lattice constant
  double rho, phi, theta; // index of the coordinate that makes a step (x,y or z)
  int M = 100000;
  int N = 100;
  int steps = 100;
  int Block_size = M / N;// throws per block
  vector <double> r(3); // coordinates vector: x y z
  // Index meaning: squaredDistancesLattice[i][m] = |r_i|^2 for the m-th lattice walk at step i.
  std::vector<std::vector<double>> squaredDistLat(steps+1, std::vector<double>(M, 0.0));
  //steps+1 the first one is zero

//---------------------main loop M random walks--------------------------
  for(int m=0; m< M; m++){
        for(int l=0 ; l<3;l++){
         r[l] = 0.0;
        }
        squaredDistLat[0][m]= 0.0;
        rho = 0;
        //performing 100 steps
        for(int i=1 ; i<steps+1;i++){
            phi = rnd.Rannyu(0., 2*M_PI);
            theta = rnd.Rannyu(0., M_PI);
            rho = rho + 2*a*(double(rnd.Head_Cross())-0.5);// step RW
            r[0] = rho*sin(theta)*cos(phi);
            r[1] = rho*sin(theta)*sin(phi);
            r[2] = rho*cos(theta);
            squaredDistLat[i][m] = get_squared_modulus(r);
        }
  }

// ---------------Blocks statistics and saving data---------------
  ofstream WriteData;
  WriteData.open("Positions_RW2_try.out");
  
  std::vector<std::vector<double>> BlocksAvg(steps+1, std::vector<double>(N, 0.0));
  std::vector<std::vector<double>> BlocksAvg2(steps+1, std::vector<double>(N, 0.0));
  vector <double> r_mean(steps+1,0.0); //mean of i step 
  vector <double> r_mean2(steps+1,0.0); //mean of i step 
  vector <double> staterror(steps+1); // error on i step
  double acc;
  for(int i=0 ; i<steps+1;i++){
    
    for(int n = 0; n < N ; n++){
      acc = 0.0;
      for(int l=0; l < Block_size; l++){
        acc += squaredDistLat[i][l+Block_size*n];
      }
      BlocksAvg[i][n] = sqrt(acc/ double(Block_size));// variable sqrt of mean of squared distance
      BlocksAvg2[i][n] = acc/ double(Block_size);
      r_mean[i] += BlocksAvg[i][n];
      r_mean2[i] += BlocksAvg2[i][n];
    }
    staterror[i] = sqrt((r_mean2[i]/double(N)-(r_mean[i] /double(N))*(r_mean[i] /double(N)))/double(N-1));
    
    WriteData << scientific << i+1  << " " <<  r_mean[i]/double(N) << " " << staterror[i] << endl;
  }
  WriteData.close();
  //saving seed
  rnd.SaveSeed();
  return 0;
}
