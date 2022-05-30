#include <cmath>

#include "compute_potentials.hpp"

double beta(const double l, const double phi){
  return pow(l,2)*phi;
}

double beta_p(const double l, const double phi){
  return pow(l,2);
}

double beta_pp(const double l, const double phi){
  return 0;
}

double solve_cubic_eqn(double a0, double a1, double a2, double init_guess){
  double guess = init_guess;
  double tol = 1e-5;
  double err = a0 + a1*guess + a2*pow(guess,2.) + pow(guess,3.);
  while(fabs(err)>tol){
    guess = guess - ((a0 + a1*guess + a2*pow(guess,2.) + pow(guess,3.))/(a1 + 2.*a2*guess + 3.*pow(guess,2.)));
    err = a0 + a1*guess + a2*pow(guess,2.) + pow(guess,3.);
  }
  return guess;

}
