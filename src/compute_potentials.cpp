#include <cmath>
#include<vector>
using std::vector;

#include "compute_potentials.hpp"
//==============================================================================
double r_of_x(double l,double x){
  return (0.5)*(l*x)/(l - x);
  // return pow(x,2)/(l*(1 - x/l));
}
double r_p_of_x(double l,double x){
  return (0.5)*pow(l,2.0)/(pow(l - x, 2.0)) ;
  // return ((2*l - x)*x)/pow(l - x,2);
}
//==============================================================================

double beta(const double l, const double phi){
  return pow(l,2)*phi;
}

double beta_p(const double l, const double phi){
  return pow(l,2);
}

double beta_pp(const double l, const double phi){
  return 0.;
}

void beta_gen(const double l, const std::vector<double> &phi_v, std::vector<double> &beta_v1,
std::vector<double> &beta_p1,
std::vector<double> &beta_pp1){

  for(int i = 0; i<phi_v.size(); i++){
    beta_v1[i] = beta(l, phi_v[i]);
    beta_p1[i] = beta_p(l, phi_v[i]);
    beta_pp1[i] = beta_pp(l, phi_v[i]);
  }
}

//==============================================================================

double solve_cubic_eqn(double a0, double a1, double a2, double a3, double init_guess){
  double guess = init_guess;
  double tol = 1e-10;
  double err = a0 + a1*guess + a2*pow(guess,2.) + a3*pow(guess,3.);
  while(fabs(err)>tol){
    guess = guess - ((a0 + a1*guess + a2*pow(guess,2.) + a3*pow(guess,3.))/(a1 + 2.*a2*guess + 3.*a3*pow(guess,2.)));
    err = a0 + a1*guess + a2*pow(guess,2.) + a3*pow(guess,3.);
  }
  return guess;

}
