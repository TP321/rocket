#ifndef RK4_H
#define RK4_H

#include <iostream>
#include <cstdlib>
#include "Parameters.h"

class RK4 {  // 4th order Runge Kutta
  private:
  double m1, v1, v2, v3, v4, h1, h2, h3, h4;   // coefficients (k1,k2,k3,k4) for RK4
  double m, v, h, t;   // variables to store mass, velocity, height and time
  double tof, tv, hmax;   // variables to store time of flight, terminal velocity and maximum altitude

  public: 
// Define velocity ODEs
  double velocity (const double g, const double rho, double v, const double Cd,
             const double A, double m, const double mf_dot, const double ue);

// Function to perform RK4 and write output on file
  void RungeKutta4 (Parameters myP_env, Parameters myP_rkt);
}; 


#endif // RK4_H
