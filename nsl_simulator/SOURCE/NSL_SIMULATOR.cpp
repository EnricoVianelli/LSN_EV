#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  string path = argv[1];
  System SYS;
  if(argc =! 1){
    cerr << "Error: missing Output dir (Where do you want the results?)." << endl;
    exit(EXIT_FAILURE);
  }
  SYS.initialize(path);
  SYS.initialize_properties();
  SYS.block_reset(0);

  for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
    for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
      SYS.step();
      SYS.measure();
      if(j%20 == 0){
        //SYS.write_XYZ(nconf); //Write actual configuration in XYZ format
        nconf++;
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();

  cout << "Results printed in "+ path << endl;

  return 0;
}
