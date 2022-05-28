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
  void solve(const Grid_data grid,Field &n_v,Field &s_v, const Field &p_v, const Field &q_v);

private:
  double rhs_shift(double r, double ss, const double P, const double Q);
  double rhs_lapse(double r, double nn, double ss, const double P, const double Q);
  void solve_shift(const Grid_data grid, Field &s_v, const Field &p_v, const Field &q_v);
  void solve_lapse(const Grid_data grid, Field &n_v,Field &s_v, const Field &p_v, const Field &q_v);
};
#endif
