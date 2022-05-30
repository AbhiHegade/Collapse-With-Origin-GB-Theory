#ifndef _COMPUTE_POTENTIALS_HPP_
#define _COMPUTE_POTENTIALS_HPP_

double beta(const double l, const double phi);

double beta_p(const double l, const double phi);

double beta_pp(const double l, const double phi);

double solve_cubic_eqn(double a0, double a1, double a2, double init_guess = 1e-5);

#endif
