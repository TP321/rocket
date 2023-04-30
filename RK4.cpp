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

using namespace std;


  double RK4::velocity (const double g, const double rho, double v, const double Cd,
             const double A, double m, const double mf_dot, const double ue) {
  double F = - g - 0.5 * rho * v * fabs(v) * (Cd * A / m) + v /fabs(v) * (mf_dot * ue / m);
  return F;
  }

// Function to perform RK4 and write output on file
  // Open file (write mode)
  void RK4::RungeKutta4(Parameters myP_env, Parameters myP_rkt){
  fstream myfile;
  myfile.open ("output.txt", ios::out);
  if (!myfile.good()){
    cout << "Error: Unable to open the output file." << endl;
  }

  // set initial conditions
  m = accumulate(myP_rkt.m0.begin(), myP_rkt.m0.end(), 0.0); 
  v = myP_env.v0;
  h = myP_env.h0;
  t = 0;
  double mf;   // fuel mass

  // Write initial conditions to the output file (space-separated)
  myfile.precision (10);
  myfile << " Time (s) " << setw(20) << " Height (m) " << setw(20) << " Velocity (m/s) " << setw(20) << " Mass (kg) " << endl; 
  myfile << t << setw(20) << h << setw(20) << v << setw(20) << m << endl; 
  
  for (int j = 0; j < myP_rkt.i; j++){  // for-loop to loop over all stages
    mf = myP_rkt.m0[j] - myP_rkt.mR[j];  // calculate fuel mass
    if (h < 0){
      break;
    }
    while (h >= 0 && mf > myP_rkt.mf_dot[j]* myP_rkt.dt[j]){   // check if there's not enough fuel for the next iteration
    
      double m1 = -myP_rkt.dt[j] * (myP_rkt.mf_dot[j]);

      double v1 = myP_rkt.dt[j] * velocity (myP_env.g, myP_env.rho, v, myP_env.Cd, myP_rkt.A[j], m, myP_rkt.mf_dot[j], myP_rkt.ue[j]);
      double v2 = myP_rkt.dt[j] * velocity (myP_env.g, myP_env.rho, v + 0.5 * v1, myP_env.Cd, myP_rkt.A[j], m + 0.5 * m1, myP_rkt.mf_dot[j], myP_rkt.ue[j]);
      double v3 = myP_rkt.dt[j] * velocity (myP_env.g, myP_env.rho, v + 0.5 * v2, myP_env.Cd, myP_rkt.A[j], m + 0.5 * m1, myP_rkt.mf_dot[j], myP_rkt.ue[j]);
      double v4 = myP_rkt.dt[j] * velocity (myP_env.g, myP_env.rho, v + v3, myP_env.Cd, myP_rkt.A[j], m + m1, myP_rkt.mf_dot[j], myP_rkt.ue[j]);

      double h1 = myP_rkt.dt[j] * (v);
      double h2 = myP_rkt.dt[j] * (v + 0.5 * v1);
      double h3 = myP_rkt.dt[j] * (v + 0.5 * v2);
      double h4 = myP_rkt.dt[j] * (v + v3);

      m += m1 / 6.0 + (m1 + m1) / 3.0 + m1 / 6.0;
      v += v1 / 6.0 + (v2 + v3) / 3.0 + v4 / 6.0;
      h += h1 / 6.0 + (h2 + h3) / 3.0 + h4 / 6.0;

  
     // update results for time and fuel mass
      t += myP_rkt.dt[j];
      mf -= myP_rkt.dt[j] * myP_rkt.mf_dot[j]; 
    
     // write values at each time step to the output file
     myfile.precision (10);
     myfile << t << setw(20) << h << setw(20) << v << setw(20) << m << endl;

    }; 
    
    if (j != myP_rkt.i-1){
      m = accumulate(myP_rkt.m0.begin()+j+1, myP_rkt.m0.end(), 0.0); // associated mass is lost when the fuel is depleted
    } 
    if (j == myP_rkt.i-1){  // compute total mass for freefall if the rocket has gone through all stages 
      m = myP_rkt.mR.back();
    } 
    if (h < 0){
      break;
    }
  myfile << t << setw(20) << h << setw(20) << v << setw(20) << m << endl;
    
  }; 

  double temp = 0; 

  //freefall
  while (h >= 0){
    double ue = 0;
    double mf_dot = 0;

    double m1 = -myP_rkt.dt.back() * (mf_dot);
  
    double v1 = myP_rkt.dt.back() * velocity (myP_env.g, myP_env.rho, v, myP_env.Cd, myP_rkt.A.back(), m, mf_dot, ue);
    double v2 = myP_rkt.dt.back() * velocity (myP_env.g, myP_env.rho, v + 0.5 * v1, myP_env.Cd, myP_rkt.A.back(), m + 0.5 * m1, mf_dot, ue);
    double v3 = myP_rkt.dt.back() * velocity (myP_env.g, myP_env.rho, v + 0.5 * v2, myP_env.Cd, myP_rkt.A.back(), m + 0.5 * m1, mf_dot, ue);
    double v4 = myP_rkt.dt.back() * velocity (myP_env.g, myP_env.rho, v + v3, myP_env.Cd, myP_rkt.A.back(), m + m1, mf_dot, ue);

    double h1 = myP_rkt.dt.back() * (v);
    double h2 = myP_rkt.dt.back() * (v + 0.5 * v1);
    double h3 = myP_rkt.dt.back() * (v + 0.5 * v2);
    double h4 = myP_rkt.dt.back() * (v + v3);

    m += m1 / 6.0 + (m1 + m1) / 3.0 + m1 / 6.0;
    v += v1 / 6.0 + (v2 + v3) / 3.0 + v4 / 6.0;
    h += h1 / 6.0 + (h2 + h3) / 3.0 + h4 / 6.0;

    t += myP_rkt.dt.back();
    myfile.precision (10);
    myfile << t << setw(20) << h << setw(20) << v << setw(20) << m << endl;

    // maximum altitude
    if (h > temp ){
      temp = h;
    }
  } 
  cout << temp << endl;

  // time of flight
  tof = t - myP_rkt.dt.back();
  cout << tof << endl;
  
  // terminal velocity
  tv = v;
  cout << tv << endl;

  // Close file
  myfile.close();

  }
 

