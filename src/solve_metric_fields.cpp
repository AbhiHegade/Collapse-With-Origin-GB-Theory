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
  double Bep, double Bepp)
  {
    double ssr = ss/r;
    double Qr = Q/r;
    double U = 0.;
    double sm1 = 0.;
    double s0 = 0.;
    double s1 = 0.;
    double s2 = 0.;
    double s3 = 0.;

    U = (-1 + 8*Bep*Qr + 8*Bep*P*ssr)*(-1 + 8*Bep*Qr + 12*Bep*P*ssr);
    sm1 = -0.25*((-1 + 8*Bep*Qr)*(pow(P,2) + pow(r,2)*pow(Qr,2)))/ssr;
    s0 = -0.5*(P*(4*Bep*pow(P,2) - pow(r,2)*Qr + 12*Bep*pow(r,2)*pow(Qr,2)));
    s1 = -4*r_Der_Q*(-Bep + 8*pow(Bep,2)*Qr)*ssr + ((-1 + 8*Bep*Qr - 10*Bep*pow(r,2)*pow(P,2)*Qr + 8*Bepp*pow(r,2)*pow(Qr,2) - 64*Bep*Bepp*pow(r,2)*pow(Qr,3) + 2*Bep*pow(r,4)*pow(Qr,3))*ssr)/2.;
    s2 = -32*pow(Bep,2)*r_Der_Q*P*pow(ssr,2) - 4*Bep*r*r_Der_P*(-1 + 8*Bep*Qr)*pow(ssr,2) - 4*(-(Bep*P) - Bepp*pow(r,2)*P*Qr + 16*Bep*Bepp*pow(r,2)*P*pow(Qr,2))*pow(ssr,2);
    s3 = -48*pow(Bep,2)*r*r_Der_P*P*pow(ssr,3) + 16*pow(Bep,2)*pow(r,2)*r_Der_Q*Qr*pow(ssr,3) + 2*Bep*r*(-(r*Qr) - 24*Bepp*r*pow(P,2)*Qr + 8*Bepp*pow(r,3)*pow(Qr,3))*pow(ssr,3);
    // if(ssr*r>1.5){
    //   cout<<U<<" "<<sm1<<" "<<s0<<" "<<s1<<" "<<s2<<" "<<s3<<endl;
    //   cout<<r_Der_Q<<" "<<Bep<<" "<<Qr<<" "<<Bepp<<" "<<P<<" "<<r_Der_P<<endl;
    //   std::exit(0);
    // }
    return (sm1/U + s0/U + s1/U + s2/U + s3/U) ;

}
//================================================================================
double Solve_metric_fields::drhs_shiftdshift(double r,
  double ss, double P, double r_Der_P, double Q, double r_Der_Q, double Bep, double Bepp)
  {

    double Qr = Q/r;
    double ssr = ss/r;
    double ans = 0.;

    ans = -2*pow(r,3)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,2) + 16*Bep*Qr*pow(r,3)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,2) - 24*Bep*Qr*pow(r,5)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,4) + pow(Qr,2)*pow(r,2)*pow(-r + 8*Bep*Qr*r,2)*(-pow(r,3) + 16*Bepp*pow(r,3)*pow(ssr,2)) + 4*Bep*pow(Qr,3)*pow(r,3)*pow(-r + 8*Bep*Qr*r,2)*(2*pow(r,2) - 32*Bepp*pow(r,2)*pow(ssr,2) + pow(r,4)*pow(ssr,2) + 48*Bepp*pow(r,4)*pow(ssr,4)) - 32*Bep*pow(r,4)*(-r + 8*Bep*Qr*r)*pow(ssr,3)*P - 32*Bepp*Qr*pow(r,6)*(-r + 8*Bep*Qr*r)*pow(ssr,3)*P - 64*pow(Bep,2)*Qr*pow(r,4)*(-r + 8*Bep*Qr*r)*pow(ssr,3)*(-4 + 5*pow(r,2)*pow(ssr,2))*P + 8*Bep*pow(Qr,2)*pow(r,4)*(-r + 8*Bep*Qr*r)*ssr*(-5*pow(r,2) + 96*Bepp*pow(r,2)*pow(ssr,2))*P + 64*pow(Bep,2)*pow(Qr,3)*pow(r,4)*(-r + 8*Bep*Qr*r)*ssr*(5*pow(r,2) - 64*Bepp*pow(r,2)*pow(ssr,2) + 40*Bepp*pow(r,4)*pow(ssr,4))*P - pow(r,5)*pow(P,2) + 24*Bep*Qr*pow(r,5)*pow(P,2) + 20*Bep*Qr*pow(r,7)*pow(ssr,2)*pow(P,2) - 128*pow(Bep,2)*pow(r,5)*pow(ssr,4)*pow(P,2) + 1024*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4)*pow(P,2) - 896*Bep*Bepp*Qr*pow(r,7)*pow(ssr,4)*pow(P,2) - 768*pow(Bep,3)*Qr*pow(r,7)*pow(ssr,6)*pow(P,2) + 192*pow(Bep,2)*pow(Qr,2)*pow(r,3)*(-pow(r,2) - 4*pow(r,4)*pow(ssr,2) + 80*Bepp*pow(r,4)*pow(ssr,4))*pow(P,2) + 128*pow(Bep,3)*pow(Qr,3)*pow(r,3)*(4*pow(r,2) + 38*pow(r,4)*pow(ssr,2) - 512*Bepp*pow(r,4)*pow(ssr,4) - 3*pow(r,6)*pow(ssr,4) + 48*Bepp*pow(r,6)*pow(ssr,6))*pow(P,2) + 40*Bep*pow(r,5)*ssr*pow(P,3) + 128*pow(Bep,2)*Qr*pow(r,3)*ssr*(-5*pow(r,2) - 3*pow(r,4)*pow(ssr,2) + 60*Bepp*pow(r,4)*pow(ssr,4))*pow(P,3) - 512*pow(Bep,3)*pow(Qr,2)*pow(r,3)*ssr*(-5*pow(r,2) - 9*pow(r,4)*pow(ssr,2) + 120*Bepp*pow(r,4)*pow(ssr,4))*pow(P,3) - 64*pow(Bep,2)*pow(r,2)*pow(ssr,2)*(7*pow(r,3) - 56*Bep*Qr*pow(r,3) - 30*Bep*Qr*pow(r,5)*pow(ssr,2) + 288*Bep*Bepp*Qr*pow(r,5)*pow(ssr,4))*pow(P,4) + 1536*pow(Bep,3)*pow(r,5)*pow(ssr,3)*pow(P,5) + r_Der_Q*(16*Bep*pow(r,3)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,2) - 128*pow(Bep,2)*Qr*pow(r,3)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,2) + 192*pow(Bep,2)*Qr*pow(r,5)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,4) + 256*pow(Bep,2)*pow(r,4)*(-r + 8*Bep*Qr*r)*pow(ssr,3)*P - 2048*pow(Bep,3)*Qr*pow(r,4)*(-r + 8*Bep*Qr*r)*pow(ssr,3)*P + 2560*pow(Bep,3)*Qr*pow(r,6)*(-r + 8*Bep*Qr*r)*pow(ssr,5)*P + 1024*pow(Bep,3)*pow(r,5)*pow(ssr,4)*pow(P,2) - 8192*pow(Bep,4)*Qr*pow(r,5)*pow(ssr,4)*pow(P,2) + 6144*pow(Bep,4)*Qr*pow(r,7)*pow(ssr,6)*pow(P,2)) + r_Der_P*(32*Bep*pow(r,4)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,3) - 256*pow(Bep,2)*Qr*pow(r,4)*pow(-r + 8*Bep*Qr*r,2)*pow(ssr,3) + 896*pow(Bep,2)*pow(r,5)*(-r + 8*Bep*Qr*r)*pow(ssr,4)*P - 7168*pow(Bep,3)*Qr*pow(r,5)*(-r + 8*Bep*Qr*r)*pow(ssr,4)*P + 7680*pow(Bep,3)*pow(r,6)*pow(ssr,5)*pow(P,2) - 61440*pow(Bep,4)*Qr*pow(r,6)*pow(ssr,5)*pow(P,2) - 18432*pow(Bep,4)*pow(r,6)*pow(ssr,6)*pow(P,3))
    ;

    ans /= 4*pow(r,2)*pow(ssr,2)*pow(pow(r,2) - 16*Bep*Qr*pow(r,2) + 64*pow(Bep,2)*pow(Qr,2)*pow(r,2) - 20*Bep*pow(r,2)*ssr*P + 160*pow(Bep,2)*Qr*pow(r,2)*ssr*P + 96*pow(Bep,2)*pow(r,2)*pow(ssr,2)*pow(P,2),2);

    return ans;
  }
//=================================================================================
double Solve_metric_fields::rhs_lapse(double r,
  double nn,
  double ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp)
  {
    if((fabs(ss)<1e-6) && ( r<10.0) ){
      return 0.;
    }
    else{
    double ssr = ss/r;
    double Qr = Q/r;
    double U = 0.;
    double sm1 = 0.;
    double s0 = 0.;
    double s1 = 0.;
    double s2 = 0.;

    U = (-1 + 8*Bep*Qr + 8*Bep*P*ssr)*(-1 + 8*Bep*Qr + 12*Bep*P*ssr);
    sm1 = (r*nn*P*Qr*(-1 + 8*Bep*Qr))/(2.*ssr);
    s0 = -(Bep*r*nn*Qr*(-5*pow(P,2) + pow(r,2)*pow(Qr,2)));
    s1 = 4*Bep*r_Der_P*nn*(-1 + 8*Bep*Qr)*ssr + 4*Bepp*r*nn*P*Qr*(-1 + 8*Bep*Qr)*ssr;
    s2 = 48*pow(Bep,2)*r_Der_P*nn*P*pow(ssr,2) - 16*pow(Bep,2)*r*r_Der_Q*nn*Qr*pow(ssr,2) - 2*Bep*nn*(-(r*Qr) - 24*Bepp*r*pow(P,2)*Qr + 8*Bepp*pow(r,3)*pow(Qr,3))*pow(ssr,2);

    return (sm1/U + s0/U + s1/U + s2/U);

  }

  }
//==============================================================================
//Interpolation and compatification functions
double Solve_metric_fields::interp_4_c(double ym2, double ym1, double y0, double y1, double y2){
  return (90*y0 + 60*y1 - 5*y2 - 20*ym1 + 3*ym2)/128.;
}
double Solve_metric_fields::interp_4_p0(double y0, double y1, double y2, double y3, double y4){
  return (35*y0 + 140*y1 - 70*y2 + 28*y3 - 5*y4)/128.;
}
double Solve_metric_fields::interp_4_pm1(double ym1, double y0, double y1, double y2, double y3){
  return (60*y0 + 90*y1 - 20*y2 + 3*y3 - 5*ym1)/128.;
}

//==================================================================================
//RK2 evolution for shift. The value at the excision point will be determined later.
//==================================================================================

void Solve_metric_fields::solve_shift(const Grid_data grid,Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v){

  double nx = grid.nx;
  int exc_i = grid.exc_i;
  vector<double> dr = grid.dr;
  vector<double> r = grid.r;
  double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;
  double dx = grid.dx;
  vector<double> x = grid.x;
  double cl = grid.cl;

  double k1 = 0.,k2=0.;
  double r_Der_P_i=0.,r_Der_Q_i = 0.,Bep_i =0.,Bepp_i = 0.;
  double r_Der_P_ip1=0.,r_Der_Q_ip1 = 0.,Bep_ip1 =0.,Bepp_ip1 = 0.;
  /*-------------------------------------------------------------------------*/

  if(exc_i==0){
    {
      int i = exc_i;
      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      if((fabs(ls)< 1e-1) || (fabs(lexp)< 1e-1)){

        if(fabs(p_v.v[i])>1e-6){
          s_v.v[i] = (fabs(p_v.v[i])/(pow(6.,0.5)))*dr[i];
        }
        else{
          s_v.v[exc_i] = 1e-10;

        }

      }

      else {
        if((fabs(p_v.v[i])>1e-6)){
        double q1 = Dx_ptpc_2nd(q_v.v[i+1], -q_v.v[i+1], dr[i]);
        double a0 = pow(p_v.v[i],2.);
        double a1 = 0.;
        double a2 = (-6. + 48.*q1*Bep_i);
        double a3 = 48.*p_v.v[i]*Bep_i;
        double root = 0;
        // root = solve_cubic_eqn((dr[i]*dr[i]*dr[i])*a0, (dr[i]*dr[i])*a1,  a2*(dr[i]), a3, (fabs(p_v.v[i])/pow(6.,0.5))*dr[i]);
        s_v.v[i] = (fabs(p_v.v[i])/(pow(6.,0.5)))*dr[i] + fabs( (2./3.)*(pow(p_v.v[i], 3))*dr[i]*Bep_i + 2.*pow((2./3.) , 0.5)*p_v.v[i]*q1*dr[i]*Bep_i );
        // s_v.v[i] = fabs(root);
      }
      else{
        s_v.v[exc_i] = 1e-10;
      }
    }
    }
    {
      int i = exc_i;
      r_Der_P_i = 0.;
      r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);

      r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], -q_v.v[i+1], dr[i]);
      r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i] + k1, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      s_v.v[i+1] = fabs(s_v.v[i] + 0.5*(k1+k2));

    }

    for(int i=exc_i+1; i<nx-2; i++){
      r_Der_P_i = Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
      r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);

      r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
      r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i] + k1, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      s_v.v[i+1] = fabs(s_v.v[i] + 0.5*(k1+k2));
    }
    s_v.v[nx-1] = 0.;
  }

  else{
    for(int i =0; i<exc_i; i++){
      s_v.v[i] = 0.;

    }
    for(int i = grid.exc_i; i< grid.exc_i+1; i++){
      r_Der_P_i = Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
      r_Der_P_ip1 = Dx_ptp0_2nd(p_v.v[i+3], p_v.v[i+2],p_v.v[i+1], dr[i+1]);

      r_Der_Q_i = Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);
      r_Der_Q_ip1 = Dx_ptp0_2nd(q_v.v[i+3], q_v.v[i+2],q_v.v[i+1], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i] + k1, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      s_v.v[i+1] = fabs(s_v.v[i] + 0.5*(k1+k2));


    }
    for(int i = grid.exc_i+1; i<nx-2; i++){
      r_Der_P_i = Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
      r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);

      r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
      r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i] + k1, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      s_v.v[i+1] = fabs(s_v.v[i] + 0.5*(k1+k2));

      // double tol = 1e-10;
      // double err = s_v.v[i+1] - s_v.v[i] - 0.5*dx*(r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i)+
      //                                             r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1));
      // while(err>tol){
      //   double fx = err;
      //   double fpx = 1. - 0.5*dx*r_p_of_x(cl, xip1)*drhs_shiftdshift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
      //
      //   s_v.v[i+1] = s_v.v[i+1] - fx/fpx;
      //   err = s_v.v[i+1] - s_v.v[i] - 0.5*dx*(r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i)+
      //                                               r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1));
      // }
    }
    // {
    //   int i = exc_i;
    //   r_Der_P_i = Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
    //   r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);
    //
    //   r_Der_Q_i = Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);
    //   r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);
    //
    //   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
    //   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
    //   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
    //   double xip1 = x[i+1];
    //   double rip1 = r_of_x(cl, xip1);
    //   k2 = dx*r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i] + k1, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
    //
    //   s_v.v[i+1] = fabs(s_v.v[i] + 0.5*(k1+k2));
    //
    //   // double tol = 1e-10;
    //   // double err = s_v.v[i+1] - s_v.v[i] - 0.5*dx*(r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i)+
    //   //                                             r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1));
    //   // while(err>tol){
    //   //   double fx = err;
    //   //   double fpx = 1. - 0.5*dx*r_p_of_x(cl, xip1)*drhs_shiftdshift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
    //   //
    //   //   s_v.v[i+1] = s_v.v[i+1] - fx/fpx;
    //   //   err = s_v.v[i+1] - s_v.v[i] - 0.5*dx*(r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i)+
    //   //                                               r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1));
    //   // }
    // }
    // for(int i = exc_i+1; i<nx-2; i++){
    //   r_Der_P_i = Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
    //   r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);
    //
    //   r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
    //   r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);
    //
    //   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
    //   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
    //   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
    //   double xip1 = x[i+1];
    //   double rip1 = r_of_x(cl, xip1);
    //   k2 = dx*r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i] + k1, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
    //
    //   s_v.v[i+1] = fabs(s_v.v[i] + 0.5*(k1+k2));
    //
    //   // double tol = 1e-10;
    //   // double err = s_v.v[i+1] - s_v.v[i] - 0.5*dx*(r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i)+
    //   //                                             r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1));
    //   // while(err>tol){
    //   //   double fx = err;
    //   //   double fpx = 1. - 0.5*dx*r_p_of_x(cl, xip1)*drhs_shiftdshift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
    //   //
    //   //   s_v.v[i+1] = s_v.v[i+1] - fx/fpx;
    //   //   err = s_v.v[i+1] - s_v.v[i] - 0.5*dx*(r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i)+
    //   //                                               r_p_of_x(cl, xip1)*rhs_shift(rip1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1));
    //   // }
    // }
    s_v.v[nx-1] = 0.;


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
  double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;
  double dx = grid.dx;
  vector<double> x = grid.x;
  double cl = grid.cl;

  double k1 = 0.,k2=0.;
  double r_Der_P_i=0.,r_Der_Q_i = 0.,Bep_i =0.,Bepp_i = 0.;
  double r_Der_P_ip1=0.,r_Der_Q_ip1 = 0.,Bep_ip1 =0.,Bepp_ip1 = 0.;
  if(exc_i==0){
    n_v.v[0] = 0.5;

    {
      int i = exc_i;
      r_Der_P_i = 0.;
      r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);

      r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], -q_v.v[i+1], dr[i]);
      r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_lapse(rip1, n_v.v[i] + k1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1 ,q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + 0.5*(k1 + k2);
    }

    for(int i =exc_i+1; i<nx-2; i++){

      r_Der_P_i = Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
      r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);

      r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
      r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_lapse(rip1, n_v.v[i] + k1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1 ,q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + 0.5*(k1 + k2);

    }

    n_v.v[nx-1] = n_v.v[nx-2];
    n_v.rescale();
  }
  else{
    for(int i=0; i<exc_i+1; i++){
      n_v.v[i] = 1;
    }
    for(int i = exc_i; i< grid.exc_i+1; i++){
      r_Der_P_i = Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
      r_Der_P_ip1 = Dx_ptp0_2nd(p_v.v[i+3], p_v.v[i+2],p_v.v[i+1], dr[i+1]);

      r_Der_Q_i = Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);
      r_Der_Q_ip1 = Dx_ptp0_2nd(q_v.v[i+3], q_v.v[i+2],q_v.v[i+1], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_lapse(rip1, n_v.v[i] + k1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1 ,q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + 0.5*(k1 + k2);

    }

    for(int i =grid.exc_i+1; i<nx-2; i++){

      r_Der_P_i = Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
      r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);

      r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
      r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);

      Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
      Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);

      Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
      Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);

      k1 = dx*r_p_of_x(cl, x[i])*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
      double xip1 = x[i+1];
      double rip1 = r_of_x(cl, xip1);
      k2 = dx*r_p_of_x(cl, xip1)*rhs_lapse(rip1, n_v.v[i] + k1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1 ,q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);

      n_v.v[i+1] = n_v.v[i] + 0.5*(k1 + k2);
      // cout<<"nn = "<<n_v.v[i+1]<<endl;

    }
    // {
    //   int i = exc_i;
    //   n_v.v[i] = 1;
    //   r_Der_P_i = Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
    //   r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);
    //
    //   r_Der_Q_i = Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);
    //   r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);
    //
    //   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
    //   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
    //   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   k1 = dx*r_p_of_x(cl, x[i])*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
    //   double xip1 = x[i+1];
    //   double rip1 = r_of_x(cl, xip1);
    //   k2 = dx*r_p_of_x(cl, xip1)*rhs_lapse(rip1, n_v.v[i] + k1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1 ,q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
    //
    //   n_v.v[i+1] = n_v.v[i] + 0.5*(k1 + k2);
    //
    // }
    //
    // for(int i =exc_i+1; i<nx-2; i++){
    //
    //   r_Der_P_i = Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
    //   r_Der_P_ip1 = Dx_ptpc_2nd(p_v.v[i+2], p_v.v[i], dr[i+1]);
    //
    //   r_Der_Q_i = Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
    //   r_Der_Q_ip1 = Dx_ptpc_2nd(q_v.v[i+2], q_v.v[i], dr[i+1]);
    //
    //   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
    //   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
    //   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
    //
    //   k1 = dx*r_p_of_x(cl, x[i])*rhs_lapse(r[i], n_v.v[i], s_v.v[i], p_v.v[i], r_Der_P_i ,q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
    //   double xip1 = x[i+1];
    //   double rip1 = r_of_x(cl, xip1);
    //   k2 = dx*r_p_of_x(cl, xip1)*rhs_lapse(rip1, n_v.v[i] + k1, s_v.v[i+1], p_v.v[i+1], r_Der_P_ip1 ,q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
    //
    //   n_v.v[i+1] = n_v.v[i] + 0.5*(k1 + k2);
    //   // cout<<"nn = "<<n_v.v[i+1]<<endl;
    //
    // }

    n_v.v[nx-1] = n_v.v[nx-2];
    n_v.rescale();
  }

}
//==============================================================================
void Solve_metric_fields::solve(const Grid_data grid,Field &n_v,Field &s_v, const Field &p_v, const Field &q_v, const Field &phi_v){
  solve_shift(grid, s_v, p_v, q_v, phi_v);
  s_v.check_isfinite(grid.t_evolve);
  solve_lapse(grid, n_v,s_v, p_v, q_v, phi_v);
  n_v.check_isfinite(grid.t_evolve);
}
//==============================================================================
//OLD Text ss RK4
// {
//   int i = exc_i;
//   r_Der_P_i = 0.;
//   r_Der_P_ip1 = Dx_ptc_4th(p_v.v[i+3], p_v.v[i+2], p_v.v[i], p_v.v[i+1], dr[i+1]);
//
//   r_Der_Q_i = Dx_ptc_4th(q_v.v[i+2], q_v.v[i+1], -q_v.v[i+1], -q_v.v[i+2], dr[i]);
//   r_Der_Q_ip1 = Dx_ptc_4th(q_v.v[i+3], q_v.v[i+2], q_v.v[i], -q_v.v[i+1], dr[i+1]);
//
//   // Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
//   // Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
//   //
//   // Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
//   // Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
//   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
//   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
//
//   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
//   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
//
//
//   // savg = (s_v.v[i] + s_v.v[i+1])/2.;
//   // pavg = (p_v.v[i] + p_v.v[i+1])/2.;
//   // qavg = (q_v.v[i] + q_v.v[i+1])/2.;
//   savg = interp_4_c(-s_v.v[i+2], -s_v.v[i+1], s_v.v[i], s_v.v[i+1], s_v.v[i+2]);
//   pavg = interp_4_c(p_v.v[i+2], p_v.v[i+1], p_v.v[i], p_v.v[i+1], p_v.v[i+2]);
//   qavg = interp_4_c(-q_v.v[i+2], -q_v.v[i+1], q_v.v[i], q_v.v[i+1], q_v.v[i+2]);
//
//   derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
//   derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;
//
//   Bepavg = (Bep_i + Bep_ip1)/2.;
//   Beppavg = (Bepp_i + Bepp_ip1)/2.;
//
//   k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
//   double xhf = x[i] + dx/2.;
//   double rhf = r_of_x(cl, xhf);
//   k2 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k1, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//
//   k3 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k2, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//   xhf = x[i+1];
//   rhf = r_of_x(cl,xhf);
//   k4 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + k3, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
//
//   s_v.v[i+1] = fabs(s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.);
//
//   if(fabs(s_v.v[i+1])<1e-2){
//     s_v.v[i+1] = fabs(s_v.v[i+1]);
//   }
//
// }
//{
// int i = exc_i+1;
//
// r_Der_P_i = Dx_ptc_4th(p_v.v[i+2], p_v.v[i+1], p_v.v[i-1], p_v.v[i], dr[i]);
// r_Der_P_ip1 = Dx_ptc_4th(p_v.v[i+3], p_v.v[i+2], p_v.v[i], p_v.v[i-1], dr[i+1]);
//
// r_Der_Q_i = Dx_ptc_4th(q_v.v[i+2], q_v.v[i+1], q_v.v[i-1], -q_v.v[i], dr[i]);
// r_Der_Q_ip1 = Dx_ptc_4th(q_v.v[i+3], q_v.v[i+2], q_v.v[i], q_v.v[i-1], dr[i+1]);
//
// Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
// Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
//
// Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
// Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
//
//
// // savg = (s_v.v[i] + s_v.v[i+1])/2.;
// // pavg = (p_v.v[i] + p_v.v[i+1])/2.;
// // qavg = (q_v.v[i] + q_v.v[i+1])/2.;
// savg = interp_4_c(-s_v.v[i+1], s_v.v[i-1], s_v.v[i], s_v.v[i+1], s_v.v[i+2]);
// pavg = interp_4_c(p_v.v[i+1], p_v.v[i-1], p_v.v[i], p_v.v[i+1], p_v.v[i+2]);
// qavg = interp_4_c(-q_v.v[i+1], q_v.v[i-1], q_v.v[i], q_v.v[i+1], q_v.v[i+2]);
//
//
// derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
// derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;
//
// Bepavg = (Bep_i + Bep_ip1)/2.;
// Beppavg = (Bepp_i + Bepp_ip1)/2.;
//
// k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
// double xhf = x[i] + dx/2.;
// double rhf = r_of_x(cl, xhf);
// k2 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k1, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//
// k3 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k2, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
// xhf = x[i+1];
// rhf = r_of_x(cl,xhf);
// k4 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + k3, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
//
//
// s_v.v[i+1] = s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
//
// if(fabs(s_v.v[i+1])<1e-2){
//   s_v.v[i+1] = fabs(s_v.v[i+1]);
// }
//}

// for(int i=exc_i+2; i<nx-3; i++){
//   r_Der_P_i = Dx_ptc_4th(p_v.v[i+2], p_v.v[i+1], p_v.v[i-1], p_v.v[i-2], dr[i]);
//   r_Der_P_ip1 = Dx_ptc_4th(p_v.v[i+3], p_v.v[i+2], p_v.v[i], p_v.v[i-1], dr[i+1]);
//
//   r_Der_Q_i = Dx_ptc_4th(q_v.v[i+2], q_v.v[i+1], q_v.v[i-1], q_v.v[i-2], dr[i]);
//   r_Der_Q_ip1 = Dx_ptc_4th(q_v.v[i+3], q_v.v[i+2], q_v.v[i], q_v.v[i-1], dr[i+1]);
//
//   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
//   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
//
//   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
//   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
//
//   // savg = (s_v.v[i] + s_v.v[i+1])/2.;
//   // pavg = (p_v.v[i] + p_v.v[i+1])/2.;
//   // qavg = (q_v.v[i] + q_v.v[i+1])/2.;
//   savg = interp_4_c(s_v.v[i-2], s_v.v[i-1], s_v.v[i], s_v.v[i+1], s_v.v[i+2]);
//   pavg = interp_4_c(p_v.v[i-2], p_v.v[i-1], p_v.v[i], p_v.v[i+1], p_v.v[i+2]);
//   qavg = interp_4_c(q_v.v[i-2], q_v.v[i-1], q_v.v[i], q_v.v[i+1], q_v.v[i+2]);
//
//   derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
//   derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;
//
//   Bepavg = (Bep_i + Bep_ip1)/2.;
//   Beppavg = (Bepp_i + Bepp_ip1)/2.;
//
//   k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
//   double xhf = x[i] + dx/2.;
//   double rhf = r_of_x(cl, xhf);
//   k2 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k1, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//
//   k3 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k2, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//   xhf = x[i+1];
//   rhf = r_of_x(cl,xhf);
//   k4 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + k3, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
//
//
//   s_v.v[i+1] = fabs(s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.);
// if(fabs(s_v.v[i+1])<1e-2){
//   s_v.v[i+1] = fabs(s_v.v[i+1]);
// }
// if(s_v.v[i+1]<-1e-20){
//   cout<<"q = "<<q_v.v[i]<<"; p = "<<p_v.v[i]<<"; phi = "<<phi_v.v[i]<<endl;
// cout<<"k1 = "<<k1<<"; k2 = "<<k2<<"; k3 = "<<k3<<"; k4 = "<<k4<<"; extra = "<<+ k1/6. + k2/3. + k3/3. + k4/6.<<endl;
// cout<<"i = "<<i<<"; s_v.v[i] = "<<s_v.v[i]<<"; s_v.v[i+1] = "<<s_v.v[i+1]<<endl;
// }

//{
//   int i = nx-3;
//
//   r_Der_P_i = Dx_ptc_4th(0., p_v.v[i+1], p_v.v[i-1], p_v.v[i-2], dr[i]);
//   r_Der_P_ip1 = Dx_ptm1_4th(0., p_v.v[i+1], p_v.v[i], p_v.v[i-1], p_v.v[i-2], dr[i+1]);
//
//   r_Der_Q_i = Dx_ptc_4th(0., q_v.v[i+1], q_v.v[i-1], q_v.v[i-2], dr[i]);
//   r_Der_Q_ip1 = Dx_ptm1_4th(0., q_v.v[i+1], q_v.v[i], q_v.v[i-1], q_v.v[i-2], dr[i+1]);
//
//   Bep_i = beta_p(ls,lexp,mu, phi_v.v[i]);
//   Bep_ip1 = beta_p(ls,lexp,mu, phi_v.v[i+1]);
//
//   Bepp_i = beta_pp(ls,lexp,mu, phi_v.v[i]);
//   Bepp_ip1 = beta_pp(ls,lexp,mu, phi_v.v[i+1]);
//
//   // savg = (s_v.v[i] + s_v.v[i+1])/2.;
//   // pavg = (p_v.v[i] + p_v.v[i+1])/2.;
//   // qavg = (q_v.v[i] + q_v.v[i+1])/2.;
//   savg = interp_4_c(s_v.v[i-2], s_v.v[i-1], s_v.v[i], s_v.v[i+1], s_v.v[i+2]);
//   pavg = interp_4_c(p_v.v[i-2], p_v.v[i-1], p_v.v[i], p_v.v[i+1], p_v.v[i+2]);
//   qavg = interp_4_c(q_v.v[i-2], q_v.v[i-1], q_v.v[i], q_v.v[i+1], q_v.v[i+2]);
//
//   derPavg = (r_Der_P_i + r_Der_P_ip1 )/2.;
//   derQavg = (r_Der_Q_i +  r_Der_Q_ip1)/2.;
//
//   Bepavg = (Bep_i + Bep_ip1)/2.;
//   Beppavg = (Bepp_i + Bepp_ip1)/2.;
//
//   k1 = dx*r_p_of_x(cl, x[i])*rhs_shift(r[i], s_v.v[i], p_v.v[i], r_Der_P_i, q_v.v[i], r_Der_Q_i, Bep_i, Bepp_i);
//   double xhf = x[i] + dx/2.;
//   double rhf = r_of_x(cl, xhf);
//   k2 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k1, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//
//   k3 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + 0.5*k2, pavg, derPavg, qavg, derQavg, Bepavg, Beppavg);
//   xhf = x[i+1];
//   rhf = r_of_x(cl,xhf);
//   k4 = dx*r_p_of_x(cl, xhf)*rhs_shift(rhf, s_v.v[i] + k3, p_v.v[i+1], r_Der_P_ip1, q_v.v[i+1], r_Der_Q_ip1, Bep_ip1, Bepp_ip1);
//
//
//   s_v.v[i+1] = s_v.v[i] + k1/6. + k2/3. + k3/3. + k4/6.;
// }

//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
