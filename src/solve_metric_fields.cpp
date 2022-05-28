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
//==============================================================================
Solve_metric_fields::Solve_metric_fields(){

}
//==============================================================================
Solve_metric_fields::~Solve_metric_fields(void)
{
}
//==============================================================================
//Note that the rhs is in d/dr. Your coordinates are in d/dx.
double Solve_metric_fields::rhs_shift(double r, double ss, const double P, const double Q){
  return
  (2*pow(r,3)*pow(P,2) + 2*pow(r,3)*pow(Q,2) + 4*pow(r,3)*P*Q*ss - 4*r*pow(ss,2))/(8.*pow(r,2)*ss)
  ;

}
//==============================================================================
double Solve_metric_fields::rhs_lapse(double r, double nn, double ss, const double P, const double Q){
  return
  -0.5*(r*nn*P*Q)/ss
  ;
}
//==================================================================================
//RK2 evolution for shift. The value at the excision point will be determined later.
//==================================================================================
void Solve_metric_fields::solve_shift(const Grid_data grid,Field &s_v, const Field &p_v, const Field &q_v){
  // assert(fabs(s_v.rbval)<1e-10);
  // double dx = grid.dx;
  double nx = grid.nx;
  int exc_i = grid.exc_i;
  vector<double> dr = grid.dr;
  vector<double> r = grid.r;

  double k1=0., k2 = 0. , k3=0., k4 = 0.;
  double  pavg = 0., qavg = 0.;
  if(exc_i==0){
    if(fabs(p_v.v[exc_i])>1e-10){
      s_v.v[exc_i] = (fabs(p_v.v[exc_i])/(pow(6.,0.5)))*r[exc_i];
    }
    else{
      s_v.v[exc_i] = 1e-10;
    }
    for(int i=exc_i; i<nx-2; i++){
      // cout<<"s_v = "<<s_v.v[i]<<endl;
      k1 = dr[i]*rhs_shift(r[i], s_v.v[i], p_v.v[i], q_v.v[i]);
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;
      k2 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + k1/2., pavg, qavg);
      k3 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + k2/2., pavg,qavg);
      k4 = dr[i]*rhs_shift(r[i] + dr[i], s_v.v[i] + k3, p_v.v[i+1], q_v.v[i+1]);
      s_v.v[i+1] = s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;

    }

    s_v.v[nx-1] = 0.;
    s_v.check_isfinite(grid.t_evolve);
  }
  else{
    for(int i =0; i<exc_i; i++){
      s_v.v[i] = 1.;

    }
    for(int i=exc_i; i<nx-2; i++){
      // cout<<"s_v = "<<s_v.v[i]<<endl;
      k1 = dr[i]*rhs_shift(r[i], s_v.v[i], p_v.v[i], q_v.v[i]);
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;
      k2 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + k1/2., pavg, qavg);
      k3 = dr[i]*rhs_shift(r[i] + 0.5*dr[i], s_v.v[i] + k2/2., pavg,qavg);
      k4 = dr[i]*rhs_shift(r[i] + dr[i], s_v.v[i] + k3, p_v.v[i+1], q_v.v[i+1]);
      s_v.v[i+1] = s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;

    }

    s_v.v[nx-1] = 0.;
    s_v.check_isfinite(grid.t_evolve);




  }


}
//==============================================================================
/*RK2 evolution for lapse. The value at the excision point is set to 1 initially
and then we rescale the solution in the end.
This works because lapse has residual gauge symmetry. */
//==================================================================================
void Solve_metric_fields::solve_lapse(const Grid_data grid,Field &n_v, Field &s_v, const Field &p_v, const Field &q_v){
  // assert(fabs(n_v.rbval-1)<1e-10);
  // double dx = grid.dx;
  double nx = grid.nx;
  int exc_i = grid.exc_i;
  vector<double> dr = grid.dr;
  vector<double> r = grid.r;

  double k1=0., k2 = 0. , k3=0., k4 = 0.;
  double nhf = 0., pavg = 0., qavg = 0., savg=0.;
  if(exc_i==0){
    n_v.v[0] = 0.5;

    for(int i =exc_i; i<nx-2; i++){
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;
      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      k1 = dr[i]*rhs_lapse(r[i], n_v.v[i], s_v.v[i],p_v.v[i], q_v.v[i]);
      k2 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k1, savg, pavg,qavg);
      k3 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k2, savg, pavg,qavg);
      k4 = dr[i]*rhs_lapse(r[i] + dr[i], n_v.v[i] + k3, s_v.v[i+1], p_v.v[i+1],q_v.v[i+1]);

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
    for(int i =exc_i; i<nx-2; i++){
      pavg = (p_v.v[i] + p_v.v[i+1])/2.;
      qavg = (q_v.v[i] + q_v.v[i+1])/2.;
      savg = (s_v.v[i] + s_v.v[i+1])/2.;
      k1 = dr[i]*rhs_lapse(r[i], n_v.v[i], s_v.v[i],p_v.v[i], q_v.v[i]);
      k2 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k1, savg, pavg,qavg);
      k3 = dr[i]*rhs_lapse(r[i] + 0.5*dr[i], n_v.v[i] + 0.5*k2, savg, pavg,qavg);
      k4 = dr[i]*rhs_lapse(r[i] + dr[i], n_v.v[i] + k3, s_v.v[i+1], p_v.v[i+1],q_v.v[i+1]);

      n_v.v[i+1] = n_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
    }
    n_v.v[nx-1] = n_v.v[nx-2];
    n_v.rescale();

    n_v.check_isfinite(grid.t_evolve);

  }

}
//==============================================================================
void Solve_metric_fields::solve(const Grid_data grid,Field &n_v,Field &s_v, const Field &p_v, const Field &q_v){
  solve_shift(grid, s_v, p_v, q_v);
  solve_lapse(grid, n_v,s_v, p_v, q_v);
}
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
