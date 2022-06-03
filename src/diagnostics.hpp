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
  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/


private:
  void find_abs_min(const vector<double> &v, double &min_elem, int &index, const double ref_val);

  int compute_radial_characteristic(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp,
  double &ingoing_c, double &outgoing_c);

  };






#endif
