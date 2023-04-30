#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
using namespace std;

class Parameters {             // class for parameters
  private:
  double rho, g, Cd, v0, h0;
  vector<double> m0;
  vector<double> mR;
  vector<double> A;
  vector<double> mf_dot;
  vector<double> ue;
  vector<double> dt;
  int i; // number of stages
    
  public:
  void getP_env (string line);

  void getP_rkt (string line);
  
  friend class RK4;
};

#endif // PARAMETERS_H
