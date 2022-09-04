#include <cmath>
#include<vector>
using std::vector;
#include<iostream>

#include "compute_potentials.hpp"
//==============================================================================
double r_of_x(double l,double x){
  return x/(1 - pow(x,2)/pow(l,2));

}
double r_p_of_x(double l,double x){
  return (pow(l,2)*(pow(l,2) + pow(x,2)))/pow(pow(l,2) - pow(x,2),2);
}
//==============================================================================

double beta(const double ls, const double lexp, const double mu, const double phi){
  double num = (1 - exp(-mu*phi*phi));
  double denom = 2.*mu;
  return (num/denom)*pow(lexp,2) + pow(ls,2)*phi ;
  // return pow(l,2)*phi;
}

double beta_p(const double ls, const double lexp, const double mu, const double phi){
  return pow(lexp,2)*phi*exp(-mu*phi*phi) + pow(ls,2);
   // return pow(l,2);
}

double beta_pp(const double ls, const double lexp, const double mu, const double phi){
  double exp1 = exp(-mu*phi*phi);
  return pow(lexp,2)*exp1*(1. - 2.*mu*phi*phi);
    // return 0.;
}

double beta_ppp(const double ls, const double lexp, const double mu, const double phi){
  double exp1 = exp(-mu*phi*phi);
  return pow(lexp,2)*exp1*(2*mu*phi*(-3 + 2*mu*pow(phi,2)));
}

void beta_gen(const double ls, const double lexp, const double mu, const std::vector<double> &phi_v, std::vector<double> &beta_v1,
std::vector<double> &beta_p1,
std::vector<double> &beta_pp1){

  for(int i = 0; i<phi_v.size(); i++){
    beta_v1[i] = beta(ls,lexp,mu,phi_v[i]);
    beta_p1[i] = beta_p(ls,lexp,mu,phi_v[i]);
    beta_pp1[i] = beta_pp(ls,lexp,mu,phi_v[i]);
  }
}

//==============================================================================

double solve_cubic_eqn(double a0, double a1, double a2, double a3, double init_guess){
  double guess = init_guess;
  double tol = 1e-10;
  int counter = 0;
  int max_steps = 10000;
  double err = a0 + a1*guess + a2*pow(guess,2.) + a3*pow(guess,3.);
  while((fabs(err)>tol)&&(counter<max_steps)){
    guess = guess - ((a0 + a1*guess + a2*pow(guess,2.) + a3*pow(guess,3.))/(a1 + 2.*a2*guess + 3.*a3*pow(guess,2.)));
    err = a0 + a1*guess + a2*pow(guess,2.) + a3*pow(guess,3.);
    counter += 1;
  }
  if(counter >= max_steps){
    std::cout<<"Max steps exceeded, stopping find root. Error = "<<err;
    std::cout<<"; polynomial coeffs, a0 = "<<a0<<"; a1 = "<<a1<<"; a2 = "<<a2<<"; a3 = "<<a3;
    std::cout<<"; init_guess = "<<init_guess<<"; root = "<<guess<<"."<<std::endl;
    if(fabs(err)>1e-6){
      std::cout<<"Error threshold exceeded, err = "<<err<<std::endl;
      std::exit(0);
    }
  }
  else{

  }
  return guess;

}
