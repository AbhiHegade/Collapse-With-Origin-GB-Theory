#ifndef _EVOLVE_SCALAR_FIELD_HPP_
#define _EVOLVE_SCALAR_FIELD_HPP_

#include <cmath>
#include <vector>
#include <string>

#include "field.hpp"
#include "grid_data.hpp"
#include "solve_metric_fields.hpp"
#include "fd_stencils.hpp"
#include "compute_potentials.hpp"

class Evolve_scalar_field{
public:
  /*-------------------Constructor--------------------------------------------*/
  Evolve_scalar_field();
  /*-------------------Destructor---------------------------------------------*/
  ~Evolve_scalar_field(void);
  //----------------------------------------------------------------------------
  void evolve(const Grid_data grid, const Field &n_v, Field &s_v, Field &p, Field &q, Field &phi);
  // void evolve_no_back_reaction(Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v);
private:
  //============================================================================
  // RHS Functions for evolution
  //============================================================================
  /* dphi/dt = rhs_phi*/
  double rhs_phi(double r,
    double nn, double r_Der_nn,
    double ss, double r_Der_ss,
    double P, double r_Der_P,
    double Q, double r_Der_Q);
  //============================================================================
  /* dq/dt = rhs_q */
  double rhs_q(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q);
  //============================================================================
  /* dp/dt = rhs_p */
  double rhs_p(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp );
  //============================================================================
  /* ds/dt = rhs_s_free, only used if excising. */
  double rhs_s_free(double r, double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp);
  //============================================================================
  /* KO filter */
  void KO_filter(const Grid_data grid,
  	const std::string type, std::vector<double> &rhs, const std::vector<double> &vec);
  //============================================================================
  void generate_rhs_non_excised(const Grid_data grid,
    const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v,
    std::vector<double> &dpdt,
    std::vector<double> &dqdt,
    std::vector<double> &dphidt);
  //============================================================================
  void generate_rhs_excised(const Grid_data grid,
    const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v,
    double &dsdt,
    std::vector<double> &dpdt,
    std::vector<double> &dqdt,
    std::vector<double> &dphidt);
  //============================================================================
  //============================================================================
  //============================================================================

};



#endif
