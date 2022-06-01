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
  /* NEED TO ADD,
  1) CHARACTERISTIC FINDER,
  2) INDEPENDENT RESIDUAL FINDER
  3) NULL - CONVERGENCE CHECK
  4) EFFECTIVE METRIC */


private:
  void find_abs_min(const vector<double> &v, double &min_elem, int &index, const double ref_val);
};
  int compute_radial_characteristic(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp,
  double &ingoing_c, double &outgoing_c);






#endif
