#include <iostream>
#include "Parameters.h"
#include "RK4.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <numeric>
#include <iomanip>
using namespace std;


int main()
{
  // read from a file to extract parameters for further analysis
  fstream myFile;
  string line1;
  string line;

  // define objects
  Parameters myP_env;
  Parameters myP_rkt;
  RK4 rk4;
    
  myFile.open ("parameters.txt", ios::in);
  // check if the file can be opened
  if (myFile.good()){ 
    while (true){  
      while(!myFile.eof()){
        getline(myFile,line1);
        myP_env.getP_env(line1); // get the enviromental parameters
          while (getline(myFile,line)){ 
            myP_rkt.getP_rkt(line);  // get the rocket_specific parameters
            }   
        }
          if (myFile.eof()){
            break;
          }
    }                      
}
  else {
        cout << "Error: Unable to open the file." << endl;
  }
  // close parameters.txt
  myFile.close();

  // perform RK4
  rk4.RungeKutta4 (myP_env, myP_rkt);

return 0;
}