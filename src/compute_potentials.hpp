#ifndef _COMPUTE_POTENTIALS_HPP_
#define _COMPUTE_POTENTIALS_HPP_

#include <vector>

double r_of_x(double cl, double x1);
double r_p_of_x(double cl, double x1);

double beta(const double ls, const double lexp, const double mu, const double phi);

double beta_p(const double ls, const double lexp, const double mu, const double phi);

double beta_pp(const double ls, const double lexp, const double mu, const double phi);

double beta_ppp(const double ls, const double lexp, const double mu, const double phi);

double solve_cubic_eqn(double a0, double a1, double a2, double a3, double init_guess = 1e-3);

void beta_gen(const double ls, const double lexp, const double mu, const std::vector<double> &phi_v, std::vector<double> &beta_v1,
std::vector<double> &beta_p1,
std::vector<double> &beta_pp1);
#endif
