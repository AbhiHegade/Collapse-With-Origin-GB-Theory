#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <cmath>
using std::pow;
using std::fabs;
#include <cassert>

#include "solve_metric_fields.hpp"
#include "field.hpp"
#include "grid_data.hpp"
#include "fd_stencils.hpp"
#include "compute_potentials.hpp"
//==============================================================================
Solve_metric_fields::Solve_metric_fields(){

}
//==============================================================================
Solve_metric_fields::~Solve_metric_fields(void)
{
}
//==============================================================================
//Note that the rhs is in d/dr. Your coordinates are in d/dx.
double Solve_metric_fields::rhs_shift(double r,
  double ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp){
  double Qr = Q/r;
  double ssr = ss/r;

  return
  (-4*Bep*r*r_Der_P*pow(ssr,2))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) + (4*r_Der_Q*ssr*(Bep - 8*pow(Bep,2)*Qr + 4*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,2) - 8*pow(Bep,2)*ssr*P))/(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - 20*Bep*ssr*P + 160*pow(Bep,2)*Qr*ssr*P + 96*pow(Bep,2)*pow(ssr,2)*pow(P,2)) + (pow(Qr,2)*pow(r,2) - 8*Bep*pow(Qr,3)*pow(r,2) - 2*pow(ssr,2) + 16*Bep*Qr*pow(ssr,2) + 16*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,2) - 128*Bep*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,2) + 4*Bep*pow(Qr,3)*pow(r,4)*pow(ssr,2) - 8*Bep*Qr*pow(r,2)*pow(ssr,4) + 64*Bep*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,4) + 2*Qr*pow(r,2)*ssr*P - 24*Bep*pow(Qr,2)*pow(r,2)*ssr*P + 16*Bep*pow(ssr,3)*P + 16*Bepp*Qr*pow(r,2)*pow(ssr,3)*P - 256*Bep*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,3)*P + pow(P,2) - 8*Bep*Qr*pow(P,2) - 20*Bep*Qr*pow(r,2)*pow(ssr,2)*pow(P,2) - 192*Bep*Bepp*Qr*pow(r,2)*pow(ssr,4)*pow(P,2) - 8*Bep*ssr*pow(P,3))/(4.*ssr*(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - 20*Bep*ssr*P + 160*pow(Bep,2)*Qr*ssr*P + 96*pow(Bep,2)*pow(ssr,2)*pow(P,2)))
  ;

}
//==============================================================================
double Solve_metric_fields::rhs_lapse(double r,
  double nn,
  double ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp)
  {
    double Qr = Q/r;
    double ssr = ss/r;

    return
    (4*Bep*r_Der_P*ssr*nn)/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) - (16*pow(Bep,2)*Qr*r*r_Der_Q*pow(ssr,2)*nn)/(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - 20*Bep*ssr*P + 160*pow(Bep,2)*Qr*ssr*P + 96*pow(Bep,2)*pow(ssr,2)*pow(P,2)) - (r*nn*(2*Bep*pow(Qr,3)*pow(r,2)*ssr - 4*Bep*Qr*pow(ssr,3) + 32*Bep*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,3) + Qr*P - 8*Bep*pow(Qr,2)*P + 8*Bepp*Qr*pow(ssr,2)*P - 64*Bep*Bepp*pow(Qr,2)*pow(ssr,2)*P - 10*Bep*Qr*ssr*pow(P,2) - 96*Bep*Bepp*Qr*pow(ssr,3)*pow(P,2)))/(2.*ssr*(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - 20*Bep*ssr*P + 160*pow(Bep,2)*Qr*ssr*P + 96*pow(Bep,2)*pow(ssr,2)*pow(P,2)))
    ;
  }

//==================================================================================
//RK2 evolution for shift. The value at the excision point will be determined later.
//==================================================================================
void Solve_metric_fields::solve_shift(const Grid_data grid,Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v){

  double nx = grid.nx;
  int exc_i = grid.exc_i;
  vector<double> dr = grid.dr;
  vector<double> r = grid.r;
  double l = grid.l;

  double k1=0., k2 = 0. , k3=0., k4 = 0.;
  double r_Der_P_i = 0., r_Der_Q_i = 0., r_Der_P_ip1 = 0., r_Der_Q_ip1 = 0.,Bep_i = 0., Bep_ip1 =0., Bepp_i = 0., Bepp_ip1 = 0.;
  double savg = 0., pavg = 0., qavg = 0.,derPavg = 0., derQavg = 0., Bepavg = 0., Beppavg=0.;

  if(exc_i==0){
    //NEED TO CHANGE THIS
    {
      int i = exc_i;
      Bep_i = beta_p(l, phi_v.v[i]);
      if(fabs(l)< 1e-5){

        if(fabs(p_v.v[i])>1e-10){
          s_v.v[i] = (fabs(p_v.v[i])/(pow(6.,0.5)))*r[i];
        }
        else{
          s_v.v[exc_i] = 1e-10;
        }

      }

      else {
        if((fabs(p_v.v[i]*Bep_i)>1e-10)){
        double q1 = (q_v.v[i+1]/dr[i]);
        double a0 = pow(p_v.v[i],2.);
        double a1 = 0.;
        double a2 = (-6. + 48.*q1*Bep_i);
        double a3 = 48.*p_v.v[i]*Bep_i;
        s_v.v[i] = fabs(solve_cubic_eqn(a0, a1, a2, a3,(fabs(p_v.v[i])/(pow(6.,0.5)))))*r[i];

      }
      else{
        s_v.v[i] = 1e-10;
      }
    }
    }
    {
      int i = exc_i;
      r_Der_P_i = 0.;
      r_Der_P_ip1 = (p_v.v[i+2] - p_v.v[i])/(2.*dr[i+1]);

      r_Der_Q_i = q_v.v[i+1]/dr[i];
      r_Der_Q_ip1 = (q_v.v[i+2] - q_v.v[i])/(2.*dr[i+1]);

      Bep_i = beta_p(l, phi_v.v[i]);
      Bep_ip1 = beta_p(l, phi_v.v[i+1]);

      Bepp_i = beta_pp(l, phi_v.v[i]);
      Bepp_ip1 = beta_pp(l, phi_v.v[i+1]);

      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;

      derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
      derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;

      Bepavg = (Bep_i + Bep_ip1)/2.;
      Beppavg = (Bepp_i + Bepp_ip1)/2.;

      k1 = dr[i]*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      k2 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + 0.5*k1, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
      k3 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + 0.5*k2, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
      k4 = dr[i]*rhs_shift(r[i] + dr[i], s_v.v[i] + k3, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      s_v.v[i+1] = s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }

    for(int i=exc_i+1; i<nx-2; i++){
      r_Der_P_i = (p_v.v[i+1] - p_v.v[i-1])/(2.*dr[i]);
      r_Der_P_ip1 = (p_v.v[i+2] - p_v.v[i])/(2.*dr[i+1]);

      r_Der_Q_i = (q_v.v[i+1] - q_v.v[i-1])/(2.*dr[i]);
      r_Der_Q_ip1 = (q_v.v[i+2] - q_v.v[i])/(2.*dr[i+1]);

      Bep_i = beta_p(l, phi_v.v[i]);
      Bep_ip1 = beta_p(l, phi_v.v[i+1]);

      Bepp_i = beta_pp(l, phi_v.v[i]);
      Bepp_ip1 = beta_pp(l, phi_v.v[i+1]);

      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;

      derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
      derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;

      Bepavg = (Bep_i + Bep_ip1)/2.;
      Beppavg = (Bepp_i + Bepp_ip1)/2.;

      k1 = dr[i]*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      k2 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + 0.5*k1, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
      k3 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + 0.5*k2, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
      k4 = dr[i]*rhs_shift(r[i] + dr[i], s_v.v[i] + k3, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      s_v.v[i+1] = s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;

    }

    s_v.v[nx-1] = 0.;
    s_v.check_isfinite(grid.t_evolve);
  }




  else{
    cout<<"Shift vector solver not implemented for excised spacetimes. Exiting."<<endl;
    std::exit(0);




  }


}
//==============================================================================
/*RK2 evolution for lapse. The value at the excision point is set to 1 initially
and then we rescale the solution in the end.
This works because lapse has residual gauge symmetry. */
//==================================================================================
void Solve_metric_fields::solve_lapse(const Grid_data grid, Field &n_v, Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v){

  double nx = grid.nx;
  int exc_i = grid.exc_i;
  vector<double> dr = grid.dr;
  vector<double> r = grid.r;
  double l = grid.l;

  double k1=0., k2 = 0. , k3=0., k4 = 0.;
  double r_Der_P_i = 0., r_Der_Q_i = 0., r_Der_P_ip1 = 0., r_Der_Q_ip1 = 0.,Bep_i = 0., Bep_ip1 =0., Bepp_i = 0., Bepp_ip1 = 0.;
  double savg = 0., pavg = 0., qavg = 0.,derPavg = 0., derQavg = 0., Bepavg = 0., Beppavg=0.;
  if(exc_i==0){
    n_v.v[0] = 0.5;

    {
      int i = exc_i;
      r_Der_P_i = 0.;
      r_Der_P_ip1 = (p_v.v[i+2] - p_v.v[i])/(2.*dr[i+1]);

      r_Der_Q_i = q_v.v[i+1]/dr[i];
      r_Der_Q_ip1 = (q_v.v[i+2] - q_v.v[i])/(2.*dr[i+1]);

      Bep_i = beta_p(l, phi_v.v[i]);
      Bep_ip1 = beta_p(l, phi_v.v[i+1]);

      Bepp_i = beta_pp(l, phi_v.v[i]);
      Bepp_ip1 = beta_pp(l, phi_v.v[i+1]);

      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;

      derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
      derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;

      Bepavg = (Bep_i + Bep_ip1)/2.;
      Beppavg = (Bepp_i + Bepp_ip1)/2.;


      k1 = dr[i]*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      k2 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k1, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k3 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k2, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k4 = dr[i]*rhs_lapse(r[i] + dr[i], n_v.v[i+1],s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1,Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }

    for(int i =exc_i+1; i<nx-2; i++){
      r_Der_P_i = (p_v.v[i+1] - p_v.v[i-1])/(2.*dr[i]);
      r_Der_P_ip1 = (p_v.v[i+2] - p_v.v[i])/(2.*dr[i+1]);

      r_Der_Q_i = (q_v.v[i+1] - q_v.v[i-1])/(2.*dr[i]);
      r_Der_Q_ip1 = (q_v.v[i+2] - q_v.v[i])/(2.*dr[i+1]);

      Bep_i = beta_p(l, phi_v.v[i]);
      Bep_ip1 = beta_p(l, phi_v.v[i+1]);

      Bepp_i = beta_pp(l, phi_v.v[i]);
      Bepp_ip1 = beta_pp(l, phi_v.v[i+1]);

      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;

      derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
      derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;

      Bepavg = (Bep_i + Bep_ip1)/2.;
      Beppavg = (Bepp_i + Bepp_ip1)/2.;


      k1 = dr[i]*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      k2 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k1, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k3 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k2, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k4 = dr[i]*rhs_lapse(r[i] + dr[i], n_v.v[i+1],s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1,Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }
    n_v.v[nx-1] = n_v.v[nx-2];
    n_v.rescale();

    n_v.check_isfinite(grid.t_evolve);
  }
  else{


    for(int i =0; i<exc_i+1; i++){
      n_v.v[i] = 0.5;

    }
    {
      int i = exc_i;
      r_Der_P_i = (-3.*p_v.v[i] + 4.*p_v.v[i+1] -p_v.v[i+2])/(2.*dr[i]);
      r_Der_P_ip1 = (p_v.v[i+2] - p_v.v[i])/(2.*dr[i+1]);

      r_Der_Q_i = (-3.*q_v.v[i] + 4.*q_v.v[i+1] -q_v.v[i+2])/(2.*dr[i]);
      r_Der_Q_ip1 = (q_v.v[i+2] - q_v.v[i])/(2.*dr[i+1]);

      Bep_i = beta_p(l, phi_v.v[i]);
      Bep_ip1 = beta_p(l, phi_v.v[i+1]);

      Bepp_i = beta_pp(l, phi_v.v[i]);
      Bepp_ip1 = beta_pp(l, phi_v.v[i+1]);

      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;

      derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
      derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;

      Bepavg = (Bep_i + Bep_ip1)/2.;
      Beppavg = (Bepp_i + Bepp_ip1)/2.;


      k1 = dr[i]*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      k2 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k1, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k3 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k2, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k4 = dr[i]*rhs_lapse(r[i] + dr[i], n_v.v[i+1],s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1,Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }

    for(int i =exc_i+1; i<nx-2; i++){
      r_Der_P_i = (p_v.v[i+1] - p_v.v[i-1])/(2.*dr[i]);
      r_Der_P_ip1 = (p_v.v[i+2] - p_v.v[i])/(2.*dr[i+1]);

      r_Der_Q_i = (q_v.v[i+1] - q_v.v[i-1])/(2.*dr[i]);
      r_Der_Q_ip1 = (q_v.v[i+2] - q_v.v[i])/(2.*dr[i+1]);

      Bep_i = beta_p(l, phi_v.v[i]);
      Bep_ip1 = beta_p(l, phi_v.v[i+1]);

      Bepp_i = beta_pp(l, phi_v.v[i]);
      Bepp_ip1 = beta_pp(l, phi_v.v[i+1]);

      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;

      derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
      derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;

      Bepavg = (Bep_i + Bep_ip1)/2.;
      Beppavg = (Bepp_i + Bepp_ip1)/2.;


      k1 = dr[i]*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      k2 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k1, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k3 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k2, savg ,pavg , derPavg, qavg, derQavg, Bepavg, Beppavg);
      k4 = dr[i]*rhs_lapse(r[i] + dr[i], n_v.v[i+1],s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1,Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }
    n_v.v[nx-1] = n_v.v[nx-2];
    n_v.rescale();

    n_v.check_isfinite(grid.t_evolve);

  }

}
//==============================================================================
void Solve_metric_fields::solve(const Grid_data grid,Field &n_v,Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v){
  solve_shift(grid, s_v, p_v, q_v, phi_v);
  solve_lapse(grid, n_v,s_v, p_v, q_v, phi_v);
}
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
