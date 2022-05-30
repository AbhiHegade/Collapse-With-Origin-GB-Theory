#ifndef _SOLVE_METRIC_FIELDS_HPP_
#define _SOLVE_METRIC_FIELDS_HPP_

#include <cmath>
#include <vector>

#include "field.hpp"
#include "grid_data.hpp"
#include "fd_stencils.hpp"

class Solve_metric_fields{
public:
  /*-------------------Constructor--------------------------------------------*/
  Solve_metric_fields();
  /*-------------------Destructor---------------------------------------------*/
  ~Solve_metric_fields(void);
  //----------------------------------------------------------------------------
  void solve(const Grid_data grid,Field &n_v,Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v);

private:
  double rhs_shift(double r,
    double ss, double P,
    double r_Der_P, double Q, double r_Der_Q, double Bep, double Bepp);
  double rhs_lapse(double r,
    double nn,
    double ss,
    double P, double r_Der_P,
    double Q, double r_Der_Q,
    double Bep, double Bepp);
  void solve_shift(const Grid_data grid, Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v);
  void solve_lapse(const Grid_data grid, Field &n_v,Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v);
};
#endif
