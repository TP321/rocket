#include "Parameters.h"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
using namespace std;

  void Parameters::getP_env (string line){  // member function to assign environment parameters
    istringstream ss (line); 
    ss >> rho >> g >> Cd >> v0 >> h0 ; 
    }

  void Parameters::getP_rkt (string line){  // member function to assign rocket-specific parameters 
    double p1, p2, p3, p4, p5, p6;
    istringstream ss (line); 
    ss >> p1 >> p2 >> p3 >> p4 >> p5 >> p6; 
    m0.push_back(p1); 
    mR.push_back(p2);
    A.push_back(p3);
    mf_dot.push_back(p4);
    ue.push_back(p5);
    dt.push_back(p6);
    i = m0.size();  // number of stages
  }
 
 

  

