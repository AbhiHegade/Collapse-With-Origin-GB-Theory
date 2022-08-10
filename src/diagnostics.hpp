#ifndef _DIAGNOSTICS_HPP_
#define _DIAGNOSTICS_HPP_

#include <cmath>
#include <vector>

#include "grid_data.hpp"
#include "field.hpp"

//==============================================================================
class Diagnostics{
public:
  /*-------------Constructor---------------------------------------------------*/
     Diagnostics();
  /*-------------Destructor----------------------------------------------------*/
     ~Diagnostics(void);
  /*---------------------------------------------------------------------------*/
  void find_apparent_horizon(Grid_data &grid, Field &s_v);
  /*---------------------------------------------------------------------------*/
  void check_for_elliptic_region(Grid_data &grid,
  const Field  &n_v, const Field &s_v,
  const Field  &p_v, const Field  &q_v,
  const Field  &phi_v,
  std::vector<double> &ingoing,
  std::vector<double> &outgoing);

  void compute_e_rr_residual(Grid_data &grid,
  const std::vector<double> &n_v,
  const std::vector<double> &s_v,
  const std::vector<double> &p_v, const std::vector<double>  &q_v,
  const std::vector<double> &phi_v,
  const std::vector<double> &n_v_np1,
  const std::vector<double> &s_v_np1, const std::vector<double> &p_v_np1,
  const std::vector<double> &q_v_np1,
  const std::vector<double> &phi_v_np1,
  std::vector<double> &residual);
  /*---------------------------------------------------------------------------*/
  void compute_NCC(Grid_data &grid,
  const std::vector<double> &n_v,
  const std::vector<double> &s_v,
  const std::vector<double> &s_v_np1,
  std::vector<double> &ncc_in,
  std::vector<double> &ncc_out);
  /*---------------------------------------------------------------------------*/
  void find_abs_min(const vector<double> &v, double &min_elem, int &index, const double ref_val, const int start_index);

  void compute_GB(Grid_data &grid,
  const std::vector<double> &n_v,
  const std::vector<double> &s_v,
  const std::vector<double> &s_v_np1,
  const std::vector<double> &n_v_np1,
  std::vector<double> &gb);

private:
  void find_outer_most_index(const vector<double> &v,int &index, const int start_index);

  int compute_radial_characteristic(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp,
  double &ingoing_c, double &outgoing_c);

  double e_rr_residual(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp,
  double t_Der_ss, double t_Der_P);

  double get_GB_Val(double r,
  double nn,
  double r_Der_nn,
  double rr_Der_nn,
  double t_Der_nn,
  double tr_Der_nn,
  double ss,
  double r_Der_ss,
  double rr_Der_ss,
  double t_Der_ss,
  double tr_Der_ss);

  double get_NCC_in(double nn,
  double r_Der_nn,
  double ss,
  double t_Der_ss );

  double get_NCC_out(double nn,
  double r_Der_nn,
  double ss,
  double t_Der_ss);

  };











#endif
