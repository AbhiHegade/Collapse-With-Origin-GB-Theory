#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <cmath>
using std::pow;
using std::fabs;
#include <cassert>
#include <string>
using std::string;

#include "field.hpp"
#include "grid_data.hpp"
#include "solve_metric_fields.hpp"
#include "evolve_scalar_field.hpp"
#include "fd_stencils.hpp"
#include "compute_potentials.hpp"
//==============================================================================
Evolve_scalar_field::Evolve_scalar_field(){

}
//==============================================================================
Evolve_scalar_field::~Evolve_scalar_field(void){

}
//==============================================================================
double Evolve_scalar_field::rhs_phi(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q){
  return nn*(P + ss*Q);
}
//==============================================================================
double Evolve_scalar_field::rhs_q(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q){

  return r_Der_P*nn + r_Der_ss*nn*Q + r_Der_Q*nn*ss + r_Der_nn*(P + Q*ss);
  // return r_Der_P*nn;
}
//==============================================================================
double Evolve_scalar_field::rhs_p(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp){

  double Qr = Q/r;
  double ssr = ss/r;
  // return
  // 2*Qr*nn + r_Der_Q*nn + r*r_Der_P*ssr*nn + r_Der_ss*nn*P + 2*ssr*nn*P + r_Der_nn*(Qr*r + r*ssr*P);
  return
  ((-64*pow(r_Der_nn,2)*(r - 8*Bep*Qr*r - 8*Bep*r*ssr*P)*(pow(Bep,2)*Qr*pow(r,5)*pow(ssr,4) + pow(Bep,2)*pow(r,3)*pow(ssr,3)*P))/nn - 64*pow(r_Der_ss,2)*(pow(Bep,2)*Qr*pow(r,4)*pow(ssr,2)*nn - 8*pow(Bep,3)*pow(Qr,2)*pow(r,4)*pow(ssr,2)*nn - 4*pow(Bep,3)*Qr*pow(r,4)*pow(ssr,3)*nn*P) - r_Der_P*(-(pow(r,5)*ssr*nn) + 16*Bep*Qr*pow(r,5)*ssr*nn - 64*pow(Bep,2)*pow(Qr,2)*pow(r,5)*ssr*nn - 32*pow(Bep,2)*pow(r,3)*pow(ssr,3)*nn + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,5)*pow(ssr,3)*nn + 32*pow(Bep,2)*pow(r,5)*pow(ssr,5)*nn - 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,7)*pow(ssr,5)*nn - 256*r_Der_Q*(-(pow(Bep,3)*pow(r,3)*pow(ssr,3)*nn) + pow(Bep,3)*pow(r,5)*pow(ssr,5)*nn) + 16*Bep*pow(r,5)*pow(ssr,2)*nn*P - 128*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,2)*nn*P - 64*pow(Bep,2)*pow(r,5)*pow(ssr,3)*nn*pow(P,2)) - r_Der_Q*(-(pow(r,4)*nn) + 16*Bep*Qr*pow(r,4)*nn - 64*pow(Bep,2)*pow(Qr,2)*pow(r,4)*nn + 16*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn + 32*pow(Bep,2)*pow(r,4)*pow(ssr,4)*nn + 16*Bep*pow(r,4)*ssr*nn*P - 96*pow(Bep,2)*Qr*pow(r,4)*ssr*nn*P + 256*pow(Bep,2)*Bepp*Qr*pow(r,4)*pow(ssr,3)*nn*P - 48*pow(Bep,2)*pow(r,4)*pow(ssr,2)*nn*pow(P,2) + 256*pow(Bep,2)*Bepp*pow(r,4)*pow(ssr,4)*nn*pow(P,2)) + (8*Qr*pow(r,4)*nn - 128*Bep*pow(Qr,2)*pow(r,4)*nn + 512*pow(Bep,2)*pow(Qr,3)*pow(r,4)*nn + 24*Bep*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn - 128*pow(Bep,2)*pow(Qr,3)*pow(r,6)*pow(ssr,2)*nn - 64*Bep*Bepp*pow(Qr,4)*pow(r,8)*pow(ssr,2)*nn + 16*Bep*pow(r,4)*pow(ssr,4)*nn - 128*Bep*Bepp*pow(Qr,2)*pow(r,6)*pow(ssr,4)*nn + 8*pow(r,4)*ssr*nn*P - 240*Bep*Qr*pow(r,4)*ssr*nn*P + 1536*pow(Bep,2)*pow(Qr,2)*pow(r,4)*ssr*nn*P - 128*Bep*Bepp*pow(Qr,3)*pow(r,6)*ssr*nn*P + 128*Bep*Bepp*Qr*pow(r,4)*pow(ssr,3)*nn*P - 128*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(ssr,3)*nn*P - 1024*Bep*pow(Bepp,2)*pow(Qr,3)*pow(r,6)*pow(ssr,3)*nn*P - 136*Bep*pow(r,4)*pow(ssr,2)*nn*pow(P,2) + 1664*pow(Bep,2)*Qr*pow(r,4)*pow(ssr,2)*nn*pow(P,2) - 64*Bep*Bepp*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn*pow(P,2) + 128*Bep*Bepp*pow(r,4)*pow(ssr,4)*nn*pow(P,2) - 1024*Bep*pow(Bepp,2)*pow(Qr,2)*pow(r,6)*pow(ssr,4)*nn*pow(P,2) + 640*pow(Bep,2)*pow(r,4)*pow(ssr,3)*nn*pow(P,3))/4. - r_Der_ss*(4*Bep*pow(Qr,2)*pow(r,6)*ssr*nn - 32*pow(Bep,2)*pow(Qr,3)*pow(r,6)*ssr*nn + 16*Bep*pow(r,4)*pow(ssr,3)*nn - 160*pow(Bep,2)*Qr*pow(r,4)*pow(ssr,3)*nn + 256*pow(Bep,2)*Bepp*pow(Qr,3)*pow(r,6)*pow(ssr,3)*nn + 256*pow(Bep,3)*Qr*pow(r,4)*r_Der_Q*pow(ssr,3)*nn - pow(r,4)*nn*P + 24*Bep*Qr*pow(r,4)*nn*P - 128*pow(Bep,2)*pow(Qr,2)*pow(r,4)*nn*P + 64*Bep*Bepp*Qr*pow(r,4)*pow(ssr,2)*nn*P - 512*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,4)*pow(ssr,2)*nn*P - 16*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn*P - 96*pow(Bep,2)*pow(r,4)*pow(ssr,4)*nn*P + 20*Bep*pow(r,4)*ssr*nn*pow(P,2) - 192*pow(Bep,2)*Qr*pow(r,4)*ssr*nn*pow(P,2) + 128*Bep*Bepp*pow(r,4)*pow(ssr,3)*nn*pow(P,2) - 1280*pow(Bep,2)*Bepp*Qr*pow(r,4)*pow(ssr,3)*nn*pow(P,2) - 80*pow(Bep,2)*pow(r,4)*pow(ssr,2)*nn*pow(P,3) - 768*pow(Bep,2)*Bepp*pow(r,4)*pow(ssr,4)*nn*pow(P,3) - 64*r_Der_P*(-(pow(Bep,2)*pow(r,3)*pow(ssr,2)*nn) + 8*pow(Bep,3)*Qr*pow(r,3)*pow(ssr,2)*nn + 2*pow(Bep,2)*pow(r,5)*pow(ssr,4)*nn - 16*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4)*nn + 4*pow(Bep,3)*pow(r,3)*pow(ssr,3)*nn*P - 12*pow(Bep,3)*pow(r,5)*pow(ssr,5)*nn*P)) - r_Der_nn*(-(Qr*pow(r,5)) + 16*Bep*pow(Qr,2)*pow(r,5) - 64*pow(Bep,2)*pow(Qr,3)*pow(r,5) - 8*Bep*pow(r,3)*pow(ssr,2) + 64*pow(Bep,2)*Qr*pow(r,3)*pow(ssr,2) + 64*Bep*Bepp*pow(Qr,2)*pow(r,5)*pow(ssr,2) - 512*pow(Bep,2)*Bepp*pow(Qr,3)*pow(r,5)*pow(ssr,2) + 4*Bep*pow(Qr,2)*pow(r,7)*pow(ssr,2) - 32*pow(Bep,2)*pow(Qr,3)*pow(r,7)*pow(ssr,2) + 16*Bep*pow(r,5)*pow(ssr,4) - 128*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,4) - pow(r,5)*ssr*P + 40*Bep*Qr*pow(r,5)*ssr*P - 256*pow(Bep,2)*pow(Qr,2)*pow(r,5)*ssr*P + 64*pow(Bep,2)*pow(r,3)*pow(ssr,3)*P + 192*Bep*Bepp*Qr*pow(r,5)*pow(ssr,3)*P - 2048*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,5)*pow(ssr,3)*P - 16*pow(Bep,2)*pow(Qr,2)*pow(r,7)*pow(ssr,3)*P - 96*pow(Bep,2)*pow(r,5)*pow(ssr,5)*P + 20*Bep*pow(r,5)*pow(ssr,2)*pow(P,2) - 256*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,2)*pow(P,2) + 128*Bep*Bepp*pow(r,5)*pow(ssr,4)*pow(P,2) - 2304*pow(Bep,2)*Bepp*Qr*pow(r,5)*pow(ssr,4)*pow(P,2) - 80*pow(Bep,2)*pow(r,5)*pow(ssr,3)*pow(P,3) - 768*pow(Bep,2)*Bepp*pow(r,5)*pow(ssr,5)*pow(P,3) - 64*pow(Bep,2)*pow(r,2)*r_Der_Q*pow(ssr,2)*(-r + 8*Bep*Qr*r + 8*Bep*r*ssr*P) - 64*r_Der_P*(-3*pow(Bep,2)*pow(r,4)*pow(ssr,3) + 24*pow(Bep,3)*Qr*pow(r,4)*pow(ssr,3) + 2*pow(Bep,2)*pow(r,6)*pow(ssr,5) - 16*pow(Bep,3)*Qr*pow(r,6)*pow(ssr,5) + 20*pow(Bep,3)*pow(r,4)*pow(ssr,4)*P - 12*pow(Bep,3)*pow(r,6)*pow(ssr,6)*P) + 16*r_Der_ss*(Bep*pow(r,3)*ssr - 16*pow(Bep,2)*Qr*pow(r,3)*ssr + 64*pow(Bep,3)*pow(Qr,2)*pow(r,3)*ssr + 8*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,3) - 64*pow(Bep,3)*pow(Qr,2)*pow(r,5)*pow(ssr,3) - 12*pow(Bep,2)*pow(r,3)*pow(ssr,2)*P + 96*pow(Bep,3)*Qr*pow(r,3)*pow(ssr,2)*P - 48*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4)*P + 32*pow(Bep,3)*pow(r,3)*pow(ssr,3)*pow(P,2))))/(pow(r,4) - 16*Bep*Qr*pow(r,4) + 64*pow(Bep,2)*pow(Qr,2)*pow(r,4) - 32*pow(Bep,2)*pow(r,4)*pow(ssr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,6)*pow(ssr,4) + 256*pow(Bep,3)*pow(r,4)*r_Der_Q*pow(ssr,4) - 16*Bep*pow(r,4)*ssr*P + 128*pow(Bep,2)*Qr*pow(r,4)*ssr*P + 64*pow(Bep,2)*pow(r,4)*pow(ssr,2)*pow(P,2) + 128*r_Der_ss*(pow(Bep,2)*pow(r,4)*pow(ssr,3) - 8*pow(Bep,3)*Qr*pow(r,4)*pow(ssr,3) - 6*pow(Bep,3)*pow(r,4)*pow(ssr,4)*P) + (128*r_Der_nn*(pow(Bep,2)*pow(r,5)*pow(ssr,4) - 8*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4) - 6*pow(Bep,3)*pow(r,5)*pow(ssr,5)*P))/nn)
  ;
}
//==============================================================================
double Evolve_scalar_field::rhs_s_free(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp){

  double Qr = Q/r;
  double ssr = ss/r;
  // return
  // (pow(Qr,2)*pow(r,3)*nn)/4. + r*r_Der_ss*ssr*nn + (r*pow(ssr,2)*nn)/2. + (Qr*r*nn*P)/(2.*ssr) + (r*nn*pow(P,2))/4.;
  return
  ((pow(Qr,2)*pow(r,3)*nn)/4. - 2*Bep*pow(Qr,3)*pow(r,3)*nn + (r*pow(ssr,2)*nn)/2. + 4*Bep*Qr*r*pow(ssr,2)*nn - 64*pow(Bep,2)*pow(Qr,2)*r*pow(ssr,2)*nn + 16*pow(Bep,2)*pow(Qr,2)*pow(r,3)*pow(ssr,4)*nn + (Qr*r*nn*P)/(2.*ssr) - (4*Bep*pow(Qr,2)*r*nn*P)/ssr + 4*Bepp*Qr*r*ssr*nn*P - 32*Bep*Bepp*pow(Qr,2)*r*ssr*nn*P - 2*Bep*pow(Qr,2)*pow(r,3)*ssr*nn*P + 4*Bep*r*pow(ssr,3)*nn*P - 128*pow(Bep,2)*Qr*r*pow(ssr,3)*nn*P + (r*nn*pow(P,2))/4. - 6*Bep*Qr*r*nn*pow(P,2) + 4*Bepp*r*pow(ssr,2)*nn*pow(P,2) - 64*Bep*Bepp*Qr*r*pow(ssr,2)*nn*pow(P,2) - 80*pow(Bep,2)*r*pow(ssr,4)*nn*pow(P,2) - 2*Bep*r*ssr*nn*pow(P,3) - 32*Bep*Bepp*r*pow(ssr,3)*nn*pow(P,3) - pow(r_Der_nn,2)*((256*pow(Bep,3)*Qr*pow(r,3)*pow(ssr,6))/nn + (256*pow(Bep,3)*r*pow(ssr,5)*P)/nn) - r_Der_P*(-4*Bep*ssr*nn + 32*pow(Bep,2)*Qr*ssr*nn + 32*pow(Bep,2)*pow(ssr,2)*nn*P) - r_Der_Q*(-4*Bep*r*pow(ssr,2)*nn + 32*pow(Bep,2)*Qr*r*pow(ssr,2)*nn + 32*pow(Bep,2)*r*pow(ssr,3)*nn*P) - pow(r_Der_ss,2)*(-128*pow(Bep,2)*r*pow(ssr,4)*nn + 768*pow(Bep,3)*Qr*r*pow(ssr,4)*nn + 768*pow(Bep,3)*r*pow(ssr,5)*nn*P) - r_Der_ss*(-(r*ssr*nn) + 12*Bep*Qr*r*ssr*nn - 32*pow(Bep,2)*pow(Qr,2)*r*ssr*nn - 16*pow(Bep,2)*pow(Qr,2)*pow(r,3)*pow(ssr,3)*nn - 256*pow(Bep,3)*r_Der_P*pow(ssr,4)*nn + 32*pow(Bep,2)*r*pow(ssr,5)*nn - 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,3)*pow(ssr,5)*nn - 256*pow(Bep,3)*r*r_Der_Q*pow(ssr,5)*nn + 12*Bep*r*pow(ssr,2)*nn*P - 96*pow(Bep,2)*Qr*r*pow(ssr,2)*nn*P - 256*pow(Bep,2)*Bepp*Qr*r*pow(ssr,4)*nn*P - 48*pow(Bep,2)*r*pow(ssr,3)*nn*pow(P,2)) - r_Der_nn*(-4*Bep*Qr*pow(r,2)*pow(ssr,2) + 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,2) - 32*pow(Bep,2)*pow(ssr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4) - 16*pow(Bep,2)*pow(Qr,2)*pow(r,4)*pow(ssr,4) + 256*pow(Bep,3)*r_Der_Q*pow(ssr,4) + 256*pow(Bep,3)*r*r_Der_P*pow(ssr,5) - 4*Bep*pow(r,2)*pow(ssr,3)*P + 32*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,3)*P + 256*pow(Bep,2)*Bepp*Qr*pow(r,2)*pow(ssr,5)*P + 16*pow(Bep,2)*pow(r,2)*pow(ssr,4)*pow(P,2) + r_Der_ss*(64*pow(Bep,2)*pow(ssr,3) - 512*pow(Bep,3)*Qr*pow(ssr,3) - 128*pow(Bep,2)*pow(r,2)*pow(ssr,5) + 1024*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,5) - 256*pow(Bep,3)*pow(ssr,4)*P + 768*pow(Bep,3)*pow(r,2)*pow(ssr,6)*P)))/(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - 32*pow(Bep,2)*pow(ssr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4) + 256*pow(Bep,3)*r_Der_Q*pow(ssr,4) - 16*Bep*ssr*P + 128*pow(Bep,2)*Qr*ssr*P + 64*pow(Bep,2)*pow(ssr,2)*pow(P,2) + r_Der_ss*(128*pow(Bep,2)*pow(ssr,3) - 1024*pow(Bep,3)*Qr*pow(ssr,3) - 768*pow(Bep,3)*pow(ssr,4)*P) + r_Der_nn*((128*pow(Bep,2)*r*pow(ssr,4))/nn - (1024*pow(Bep,3)*Qr*r*pow(ssr,4))/nn - (768*pow(Bep,3)*r*pow(ssr,5)*P)/nn))
  ;
}
//==============================================================================
void Evolve_scalar_field::KO_filter(Grid_data grid,
   const string type, vector<double> &rhs, const vector<double> &vec)
{ //Fourth Order Filter
   double eps= 0.5;
   double dt = grid.dt;
   int exc_i = grid.exc_i;
   int nx = grid.nx;

   for (int i=exc_i+3; i<nx-3; ++i) {
      rhs[i]+= (eps/(64.0*dt))*(
	 vec[i+3]
      -	 6*vec[i+2]
      +	 15*vec[i+1]
      -	 20*vec[i]
      +	 15*vec[i-1]
      -	 6*vec[i-2]
      +	 vec[i-3]
      );
   }
/*---------------------------------------------------------------------------*/
   if (exc_i>0) return;
/*---------------------------------------------------------------------------*/
   if (type=="even") {
      rhs[2]+= (eps/(64.0*dt))*(
	 vec[5]
      -	 6*vec[4]
      +	 15*vec[3]
      -	 20*vec[2]
      +	 15*vec[1]
      -	 6*vec[0]
      +	 vec[1]
      );
      rhs[1]+= (eps/(64*dt))*(
         vec[4]
      -  6*vec[3]
      +  15*vec[2]
      -  20*vec[1]
      +  15*vec[0]
      -  6*vec[1]
      +  vec[2]
      );
      rhs[0]+= (eps/(64*dt))*(
         vec[3]
      -  6*vec[2]
      +  15*vec[1]
      -  20*vec[0]
      +  15*vec[1]
      -  6*vec[2]
      +  vec[3]
      );
   } else
   if (type=="odd") {
      rhs[2]+= (eps/(64*dt))*(
         vec[5]
      -	 6*vec[4]
      +	 15*vec[3]
      -	 20*vec[2]
      +	 15*vec[1]
      -	 6*vec[0]
      +	 (-vec[1])
      );
      rhs[1]+= (eps/(64*dt))*(
         vec[4]
      -  6*vec[3]
      +  15*vec[2]
      -  20*vec[1]
      +  15*vec[0]
      -  (-6*vec[1])
      +  (-vec[2])
      );
      rhs[0]+= (eps/(64*dt))*(
         vec[3]
      -	 6*vec[2]
      +	 15*vec[1]
      -	 20*vec[0]
      +	 (-15*vec[1])
      -	 (-6*vec[2])
      +	 (-vec[3])
      );
   } else {
      /* do nothing */
   }
}
//==============================================================================
//==============================================================================
void Evolve_scalar_field::generate_rhs_non_excised(Grid_data grid,
  const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v,
  vector<double> &dpdt,
  vector<double> &dqdt,
  vector<double> &dphidt){

  //Field ordering p_v,q_v,phi_v
  assert(grid.exc_i ==0);
  vector<double> r = grid.r;
  int nx = grid.nx;
  vector<double> dr = grid.dr;
  double l = grid.l;

  double r_Der_nn=0., r_Der_ss = 0., r_Der_P = 0., r_Der_Q = 0.;
  double Bep = 0., Bepp = 0.;

  {
      int i =0;
      double nn = 0.,rr_Der_nn = 0., ss=0., P = 0., rr_Der_P = 0, Q=0.;

      nn = n_v.v[i]; r_Der_nn = 0.; rr_Der_nn = (2.*n_v.v[i+1] - 2.*n_v.v[i])/(pow(dr[i],2.));
      ss = s_v.v[i]; r_Der_ss = s_v.v[i+1]/dr[i];

      P = p_v.v[i]; r_Der_P = 0.; rr_Der_P = (2.*p_v.v[i+1] - 2.*p_v.v[i])/(pow(dr[i],2.));

      Q = q_v.v[i]; r_Der_Q = q_v.v[i+1]/dr[i];


      Bep = beta_p(l, phi_v.v[i]);
      Bepp = beta_pp(l, phi_v.v[i]);

      dpdt[i] = (-8*Bep*rr_Der_nn*pow(r_Der_ss,2)*nn + 128*pow(Bep,2)*rr_Der_nn*r_Der_Q*pow(r_Der_ss,2)*nn - 512*pow(Bep,3)*rr_Der_nn*pow(r_Der_Q,2)*pow(r_Der_ss,2)*nn + 3*r_Der_Q*pow(nn,2) - 48*Bep*pow(r_Der_Q,2)*pow(nn,2) + 192*pow(Bep,2)*pow(r_Der_Q,3)*pow(nn,2) - 32*pow(Bep,2)*rr_Der_P*pow(r_Der_ss,3)*pow(nn,2) + 256*pow(Bep,3)*rr_Der_P*r_Der_Q*pow(r_Der_ss,3)*pow(nn,2) - 12*Bep*pow(r_Der_ss,4)*pow(nn,2) + 64*pow(Bep,2)*r_Der_Q*pow(r_Der_ss,4)*pow(nn,2) + 256*pow(Bep,3)*pow(r_Der_Q,2)*pow(r_Der_ss,4)*pow(nn,2) + 128*pow(Bep,2)*rr_Der_nn*pow(r_Der_ss,3)*nn*P - 1024*pow(Bep,3)*rr_Der_nn*r_Der_Q*pow(r_Der_ss,3)*nn*P + 3*r_Der_ss*pow(nn,2)*P - 100*Bep*r_Der_Q*r_Der_ss*pow(nn,2)*P + 608*pow(Bep,2)*pow(r_Der_Q,2)*r_Der_ss*pow(nn,2)*P - 32*Bep*Bepp*r_Der_Q*pow(r_Der_ss,3)*pow(nn,2)*P + 256*pow(Bep,2)*Bepp*pow(r_Der_Q,2)*pow(r_Der_ss,3)*pow(nn,2)*P + 256*pow(Bep,3)*rr_Der_P*pow(r_Der_ss,4)*pow(nn,2)*P + 96*pow(Bep,2)*pow(r_Der_ss,5)*pow(nn,2)*P + 256*pow(Bep,3)*r_Der_Q*pow(r_Der_ss,5)*pow(nn,2)*P - 512*pow(Bep,3)*rr_Der_nn*pow(r_Der_ss,4)*nn*pow(P,2) - 54*Bep*pow(r_Der_ss,2)*pow(nn,2)*pow(P,2) + 656*pow(Bep,2)*r_Der_Q*pow(r_Der_ss,2)*pow(nn,2)*pow(P,2) - 96*Bep*Bepp*pow(r_Der_ss,4)*pow(nn,2)*pow(P,2) + 1024*pow(Bep,2)*Bepp*r_Der_Q*pow(r_Der_ss,4)*pow(nn,2)*pow(P,2) + 240*pow(Bep,2)*pow(r_Der_ss,3)*pow(nn,2)*pow(P,3) + 768*pow(Bep,2)*Bepp*pow(r_Der_ss,5)*pow(nn,2)*pow(P,3))/(nn*(1 - 16*Bep*r_Der_Q + 64*pow(Bep,2)*pow(r_Der_Q,2) + 96*pow(Bep,2)*pow(r_Der_ss,4) - 768*pow(Bep,3)*r_Der_Q*pow(r_Der_ss,4) - 16*Bep*r_Der_ss*P + 128*pow(Bep,2)*r_Der_Q*r_Der_ss*P - 768*pow(Bep,3)*pow(r_Der_ss,5)*P + 64*pow(Bep,2)*pow(r_Der_ss,2)*pow(P,2)));

      dqdt[i] = 0.;

      dphidt[i] = n_v.v[i]*p_v.v[i];
     }


  for(int i = 1; i<nx-1; i++){
    r_Der_nn = (n_v.v[i+1] - n_v.v[i-1])/(2*dr[i]);
    r_Der_ss = (s_v.v[i+1] - s_v.v[i-1])/(2*dr[i]);
    r_Der_P = (p_v.v[i+1] - p_v.v[i-1])/(2*dr[i]);
    r_Der_Q = (q_v.v[i+1] - q_v.v[i-1])/(2*dr[i]);
    Bep = beta_p(l,phi_v.v[i]);
    Bepp = beta_pp(l,phi_v.v[i]);

      dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q, Bep, Bepp);
      dqdt[i] = rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
      dphidt[i] = rhs_phi(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);


  }

  {
    int i = nx-1;
    dpdt[i] = 0.;
    dqdt[i] = 0.;
    dphidt[i] = 0.;

  }

  KO_filter(grid, "even", dpdt, p_v.v);
  KO_filter(grid, "odd", dqdt, q_v.v);
  KO_filter(grid, "even" ,dphidt, phi_v.v);

}
//==============================================================================
// Generate RHS With Excision
//==============================================================================
// void Evolve_scalar_field::generate_rhs_excised(const Grid_data grid,
//   const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v,
//   double &dsdt,
//   vector<double> &dpdt,
//   vector<double> &dqdt,
//   vector<double> &dphidt){
//
//   //Field ordering s_v, p_v,q_v,phi_v
//   assert(grid.exc_i >0);
//   vector<double> r = grid.r;
//   int nx = grid.nx;
//   vector<double> dr = grid.dr;
//
//   double r_Der_nn=0., r_Der_ss = 0., r_Der_P = 0., r_Der_Q = 0.;
//
//
//   {
//     int i = grid.exc_i;
//     r_Der_nn = (-3.*n_v.v[i] + 4.*n_v.v[i+1] - n_v.v[i+2])/(2*dr[i]);
//     r_Der_ss = (-3.*s_v.v[i] + 4.*s_v.v[i+1] - s_v.v[i+2])/(2*dr[i]);
//     r_Der_P =  (-3.*p_v.v[i] + 4.*p_v.v[i+1] - p_v.v[i+2])/(2*dr[i]);;
//     r_Der_Q =  (-3.*q_v.v[i] + 4.*q_v.v[i+1] - q_v.v[i+2])/(2*dr[i]);;
//
//     dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//     dqdt[i] = rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//     dphidt[i] = rhs_phi(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//     dsdt = rhs_s_free(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//
//
//   }
//
//
//   for(int i = grid.exc_i + 1; i<nx-1; i++){
//     r_Der_nn = (n_v.v[i+1] - n_v.v[i-1])/(2*dr[i]);
//     r_Der_ss = (s_v.v[i+1] - s_v.v[i-1])/(2*dr[i]);
//     r_Der_P = (p_v.v[i+1] - p_v.v[i-1])/(2*dr[i]);
//     r_Der_Q = (q_v.v[i+1] - q_v.v[i-1])/(2*dr[i]);
//
//       dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//       dqdt[i] = rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//       dphidt[i] = rhs_phi(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//
//
//   }
//
//   {
//     int i = nx-1;
//     dpdt[i] = 0.;
//     dqdt[i] = 0.;
//     dphidt[i] = 0.;
//
//   }
//
//   KO_filter(grid, "even", dpdt, p_v.v);
//   KO_filter(grid, "odd", dqdt, q_v.v);
//   KO_filter(grid, "even" ,dphidt, phi_v.v);
//
// }

//==============================================================================
//==============================================================================
void Evolve_scalar_field::evolve(const Grid_data grid, const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
  Solve_metric_fields solve_metric_fields;

  if(grid.exc_i>0){
    cout<<"Excision not implemented for evolution equations. Exit. "<<endl;
    std::exit(0);
    // vector<double> r = grid.r;
    // int nx = grid.nx;
    // double dt = grid.dt;
    // vector<double> dr = grid.dr;
    // int exc_i = grid.exc_i;
    //
    // vector<double> kp1(nx,0), kq1(nx,0), kphi1(nx,0);
    // vector<double> kp2(nx,0), kq2(nx,0), kphi2(nx,0);
    // vector<double> kp3(nx,0), kq3(nx,0), kphi3(nx,0);
    // vector<double> kp4(nx,0), kq4(nx,0), kphi4(nx,0);
    // double sk1 = 0, sk2 = 0, sk3 = 0, sk4 = 0;
    // vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);
    // double dsdt;
    //
    // generate_rhs_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dsdt, dpdt, dqdt, dphidt);
    //
    // for(int i =exc_i; i<nx-1; i++){
    //   kp1[i] = dt*dpdt[i];
    //   kq1[i] = dt*dqdt[i];
    //   kphi1[i] = dt*dphidt[i];
    //
    // }
    // sk1 = dt*dsdt;
    // {
    //   Field n_1("n_k1", "even", grid);
    //   Field s_1("s_k1", "odd", grid);
    //   Field p_1("p_k1", "even", grid);
    //   Field q_1("q_k1", "odd", grid);
    //   Field phi_1("phi_k1", "even", grid);
    //
    //   for(int i = exc_i; i<nx-1; i++){
    //     n_1.v[i] = n_v.v[i];
    //     s_1.v[i] = s_v.v[i] ;
    //     p_1.v[i] = p_v.v[i] + 0.5*kp1[i];
    //     q_1.v[i] = q_v.v[i] + 0.5*kq1[i];
    //     phi_1.v[i] = phi_v.v[i] + 0.5*kphi1[i];
    //   }
    //   s_1.v[exc_i] = s_v.v[exc_i] + 0.5*sk1;
    //   solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1);
    //   generate_rhs_excised(grid, n_1, s_1, p_1, q_1, phi_1, dsdt, dpdt, dqdt, dphidt);
    //
    //   for(int i =0; i<nx-1; i++){
    //     kp2[i] = dt*dpdt[i];
    //     kq2[i] = dt*dqdt[i];
    //     kphi2[i] = dt*dphidt[i];
    //
    //   }
    //   sk2 = dt*dsdt;
    // }
    //
    //   {
    //     Field n_2("n_k2", "even", grid);
    //     Field s_2("s_k2", "odd", grid);
    //     Field p_2("p_k2", "even", grid);
    //     Field q_2("q_k2", "odd", grid);
    //     Field phi_2("phi_k2", "even", grid);
    //
    //     for(int i = exc_i; i<nx-1; i++){
    //       n_2.v[i] = n_v.v[i];
    //       s_2.v[i] = s_v.v[i] ;
    //       p_2.v[i] = p_v.v[i] + 0.5*kp2[i];
    //       q_2.v[i] = q_v.v[i] + 0.5*kq2[i];
    //       phi_2.v[i] = phi_v.v[i] + 0.5*kphi2[i];
    //     }
    //     s_2.v[exc_i] = s_v.v[exc_i] + 0.5*sk2;
    //     solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2);
    //     generate_rhs_excised(grid, n_2, s_2, p_2, q_2, phi_2, dsdt, dpdt, dqdt, dphidt);
    //
    //     for(int i =exc_i; i<nx-1; i++){
    //       kp3[i] = dt*dpdt[i];
    //       kq3[i] = dt*dqdt[i];
    //       kphi3[i] = dt*dphidt[i];
    //
    //     }
    //     sk3 = dt*dsdt;
    //
    //
    //
    //   }
    //   {
    //     Field n_3("n_k3", "even", grid);
    //     Field s_3("s_k3", "odd", grid);
    //     Field p_3("p_k3", "even", grid);
    //     Field q_3("q_k3", "odd", grid);
    //     Field phi_3("phi_k3", "even", grid);
    //
    //     for(int i = exc_i; i<nx-1; i++){
    //       n_3.v[i] = n_v.v[i];
    //       s_3.v[i] = s_v.v[i];
    //       p_3.v[i] = p_v.v[i] + kp3[i];
    //       q_3.v[i] = q_v.v[i] + kq3[i];
    //       phi_3.v[i] = phi_v.v[i] + kphi3[i];
    //     }
    //     s_3.v[exc_i] = s_v.v[exc_i] + sk3;
    //     solve_metric_fields.solve(grid, n_3, s_3, p_3, q_3);
    //     generate_rhs_excised(grid, n_3, s_3, p_3, q_3, phi_3, dsdt, dpdt, dqdt, dphidt);
    //
    //     for(int i =exc_i; i<nx-1; i++){
    //       kp4[i] = dt*dpdt[i];
    //       kq4[i] = dt*dqdt[i];
    //       kphi4[i] = dt*dphidt[i];
    //
    //     }
    //     sk4 = dt*dsdt;
    //
    //
    //
    //   }
    //   for(int i=exc_i; i<nx-1; i++){
    //     p_v.v[i] += (1./6.)*kp1[i] + (1./3.)*kp2[i] + (1./3.)*kp3[i] + (1./6.)*kp4[i];
    //     q_v.v[i] += (1./6.)*kq1[i] + (1./3.)*kq2[i] + (1./3.)*kq3[i] + (1./6.)*kq4[i];
    //     phi_v.v[i] += (1./6.)*kphi1[i] + (1./3.)*kphi2[i] + (1./3.)*kphi3[i] + (1./6.)*kphi4[i];
    //   }
    //   for(int i=0; i<exc_i; i++){
    //     p_v.v[i] = p_v.v[exc_i];
    //     q_v.v[i] = s_v.v[exc_i];
    //     phi_v.v[i] = phi_v.v[exc_i];
    //   }
    //   s_v.v[exc_i] += (1./6.)*sk1 + (1./3.)*sk2 + (1./3.)*sk3 + (1./6.)*sk4;
    //
    //   p_v.check_isfinite(grid.t_evolve);
    //   q_v.check_isfinite(grid.t_evolve);
    //   phi_v.check_isfinite(grid.t_evolve);

  }
  else{
    //RK4
    assert(grid.exc_i ==0);
    vector<double> r = grid.r;
    int nx = grid.nx;
    double dt = grid.dt;
    vector<double> dr = grid.dr;

    vector<double> kp1(nx,0), kq1(nx,0), kphi1(nx,0);
    vector<double> kp2(nx,0), kq2(nx,0), kphi2(nx,0);
    vector<double> kp3(nx,0), kq3(nx,0), kphi3(nx,0);
    vector<double> kp4(nx,0), kq4(nx,0), kphi4(nx,0);
    vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);

    generate_rhs_non_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dpdt, dqdt, dphidt);

    for(int i =0; i<nx-1; i++){
      kp1[i] = dt*dpdt[i];
      kq1[i] = dt*dqdt[i];
      kphi1[i] = dt*dphidt[i];

    }
    {
      Field n_1("n_k1", "even", grid);
      Field s_1("s_k1", "odd", grid);
      Field p_1("p_k1", "even", grid);
      Field q_1("q_k1", "odd", grid);
      Field phi_1("phi_k1", "even", grid);

      for(int i = 0; i<nx-1; i++){
        n_1.v[i] = n_v.v[i];
        s_1.v[i] = s_v.v[i];
        p_1.v[i] = p_v.v[i] + 0.5*kp1[i];
        q_1.v[i] = q_v.v[i] + 0.5*kq1[i];
        phi_1.v[i] = phi_v.v[i] + 0.5*kphi1[i];
      }
      solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1, phi_1);
      generate_rhs_non_excised(grid, n_1, s_1, p_1, q_1, phi_1, dpdt, dqdt, dphidt);

      for(int i =0; i<nx-1; i++){
        kp2[i] = dt*dpdt[i];
        kq2[i] = dt*dqdt[i];
        kphi2[i] = dt*dphidt[i];

      }
    }

      {
        Field n_2("n_k2", "even", grid);
        Field s_2("s_k2", "odd", grid);
        Field p_2("p_k2", "even", grid);
        Field q_2("q_k2", "odd", grid);
        Field phi_2("phi_k2", "even", grid);

        for(int i = 0; i<nx-1; i++){
          n_2.v[i] = n_v.v[i];
          s_2.v[i] = s_v.v[i];
          p_2.v[i] = p_v.v[i] + 0.5*kp2[i];
          q_2.v[i] = q_v.v[i] + 0.5*kq2[i];
          phi_2.v[i] = phi_v.v[i] + 0.5*kphi2[i];
        }
        solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2, phi_2);
        generate_rhs_non_excised(grid, n_2, s_2, p_2, q_2, phi_2, dpdt, dqdt, dphidt);

        for(int i =0; i<nx-1; i++){
          kp3[i] = dt*dpdt[i];
          kq3[i] = dt*dqdt[i];
          kphi3[i] = dt*dphidt[i];

        }



      }
      {
        Field n_3("n_k3", "even", grid);
        Field s_3("s_k3", "odd", grid);
        Field p_3("p_k3", "even", grid);
        Field q_3("q_k3", "odd", grid);
        Field phi_3("phi_k3", "even", grid);

        for(int i = 0; i<nx-1; i++){
          n_3.v[i] = n_v.v[i];
          s_3.v[i] = s_v.v[i];
          p_3.v[i] = p_v.v[i] + kp3[i];
          q_3.v[i] = q_v.v[i] + kq3[i];
          phi_3.v[i] = phi_v.v[i] + kphi3[i];
        }
        solve_metric_fields.solve(grid, n_3, s_3, p_3, q_3, phi_3);
        generate_rhs_non_excised(grid, n_3, s_3, p_3, q_3, phi_3, dpdt, dqdt, dphidt);

        for(int i =0; i<nx-1; i++){
          kp4[i] = dt*dpdt[i];
          kq4[i] = dt*dqdt[i];
          kphi4[i] = dt*dphidt[i];

        }



      }
      for(int i=0; i<nx-1; i++){
        p_v.v[i] += (1./6.)*kp1[i] + (1./3.)*kp2[i] + (1./3.)*kp3[i] + (1./6.)*kp4[i];
        q_v.v[i] += (1./6.)*kq1[i] + (1./3.)*kq2[i] + (1./3.)*kq3[i] + (1./6.)*kq4[i];
        phi_v.v[i] += (1./6.)*kphi1[i] + (1./3.)*kphi2[i] + (1./3.)*kphi3[i] + (1./6.)*kphi4[i];
      }

      p_v.check_isfinite(grid.t_evolve);
      q_v.check_isfinite(grid.t_evolve);
      phi_v.check_isfinite(grid.t_evolve);




  }
}
//==============================================================================
// void Evolve_scalar_field::evolve_no_back_reaction(Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
//   // Solve_metric_fields solve_metric_fields(grid);
//
//   if(grid.exc_i>0){
//     cout<<"not implemented"<<endl;
//     std::exit(0);
//   }
//   else{
//     //RK4
//     assert(grid.exc_i ==0);
//     vector<double> r = grid.r;
//     int nx = grid.nx;
//     double dt = grid.dt;
//     vector<double> dr = grid.dr;
//
//     vector<double> kp1(nx,0), kq1(nx,0), kphi1(nx,0);
//     vector<double> kp2(nx,0), kq2(nx,0), kphi2(nx,0);
//     vector<double> kp3(nx,0), kq3(nx,0), kphi3(nx,0);
//     vector<double> kp4(nx,0), kq4(nx,0), kphi4(nx,0);
//     vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);
//
//     generate_rhs_non_excised(n_v, s_v, p_v, q_v, phi_v, dpdt, dqdt, dphidt);
//
//     for(int i =0; i<nx-1; i++){
//       kp1[i] = dt*dpdt[i];
//       kq1[i] = dt*dqdt[i];
//       kphi1[i] = dt*dphidt[i];
//
//     }
//     {
//       Field n_1("n_k1", "even", grid);
//       Field s_1("s_k1", "odd", grid);
//       Field p_1("p_k1", "even", grid);
//       Field q_1("q_k1", "odd", grid);
//       Field phi_1("phi_k1", "even", grid);
//
//       for(int i = 0; i<nx-1; i++){
//         n_1.v[i] = n_v.v[i];
//         s_1.v[i] = s_v.v[i];
//         p_1.v[i] = p_v.v[i] + 0.5*kp1[i];
//         q_1.v[i] = q_v.v[i] + 0.5*kq1[i];
//         phi_1.v[i] = phi_v.v[i] + 0.5*kphi1[i];
//       }
//       // solve_metric_fields.solve(n_1, s_1, p_1, q_1);
//       generate_rhs_non_excised(n_1, s_1, p_1, q_1, phi_1, dpdt, dqdt, dphidt);
//
//       for(int i =0; i<nx-1; i++){
//         kp2[i] = dt*dpdt[i];
//         kq2[i] = dt*dqdt[i];
//         kphi2[i] = dt*dphidt[i];
//
//       }
//     }
//
//       {
//         Field n_2("n_k2", "even", grid);
//         Field s_2("s_k2", "odd", grid);
//         Field p_2("p_k2", "even", grid);
//         Field q_2("q_k2", "odd", grid);
//         Field phi_2("phi_k2", "even", grid);
//
//         for(int i = 0; i<nx-1; i++){
//           n_2.v[i] = n_v.v[i];
//           s_2.v[i] = s_v.v[i];
//           p_2.v[i] = p_v.v[i] + 0.5*kp2[i];
//           q_2.v[i] = q_v.v[i] + 0.5*kq2[i];
//           phi_2.v[i] = phi_v.v[i] + 0.5*kphi2[i];
//         }
//         // solve_metric_fields.solve(n_2, s_2, p_2, q_2);
//         generate_rhs_non_excised(n_2, s_2, p_2, q_2, phi_2, dpdt, dqdt, dphidt);
//
//         for(int i =0; i<nx-1; i++){
//           kp3[i] = dt*dpdt[i];
//           kq3[i] = dt*dqdt[i];
//           kphi3[i] = dt*dphidt[i];
//
//         }
//
//
//
//       }
//       {
//         Field n_3("n_k3", "even", grid);
//         Field s_3("s_k3", "odd", grid);
//         Field p_3("p_k3", "even", grid);
//         Field q_3("q_k3", "odd", grid);
//         Field phi_3("phi_k3", "even", grid);
//
//         for(int i = 0; i<nx-1; i++){
//           n_3.v[i] = n_v.v[i];
//           s_3.v[i] = s_v.v[i];
//           p_3.v[i] = p_v.v[i] + kp3[i];
//           q_3.v[i] = q_v.v[i] + kq3[i];
//           phi_3.v[i] = phi_v.v[i] + kphi3[i];
//         }
//         // solve_metric_fields.solve(n_3, s_3, p_3, q_3);
//         generate_rhs_non_excised(n_3, s_3, p_3, q_3, phi_3, dpdt, dqdt, dphidt);
//
//         for(int i =0; i<nx-1; i++){
//           kp4[i] = dt*dpdt[i];
//           kq4[i] = dt*dqdt[i];
//           kphi4[i] = dt*dphidt[i];
//
//         }
//
//
//
//       }
//       for(int i=0; i<nx-1; i++){
//         p_v.v[i] = p_v.v[i] + (1./6.)*kp1[i] + (1./3.)*kp2[i] + (1./3.)*kp3[i] + (1./6.)*kp4[i];
//         q_v.v[i] = q_v.v[i] + (1./6.)*kq1[i] + (1./3.)*kq2[i] + (1./3.)*kq3[i] + (1./6.)*kq4[i];
//         phi_v.v[i] = phi_v.v[i] +  (1./6.)*kphi1[i] + (1./3.)*kphi2[i] + (1./3.)*kphi3[i] + (1./6.)*kphi4[i];
//       }
//
//       p_v.check_isfinite(grid.t_evolve);
//       q_v.check_isfinite(grid.t_evolve);
//       phi_v.check_isfinite(grid.t_evolve);
//
//
//
//
//   }
// }
// void Evolve_scalar_field::evolve_no_back_reaction( Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
//   // Solve_metric_fields solve_metric_fields(grid);
//
//   if(grid.exc_i>0){
//     cout<<"not implemented"<<endl;
//     std::exit(0);
//   }
//   else{
//     double exc_i = grid.exc_i;
//     assert(exc_i==0);
//     vector<double> r = grid.r;
//     int nx = grid.nx;
//     double dt = grid.dt;
//     vector<double> dr = grid.dr;
//     vector<double> p_np1(nx,0), q_np1(nx,0), phi_np1(nx,0);
//     /*========================================================================*/
//     //RK1 for testing.
//     /*========================================================================*/
//      double r_Der_nn=0., r_Der_ss = 0., r_Der_P = 0., r_Der_Q = 0.;
//      n_v.v[0] = (4./3.)*n_v.v[1] - (1./3.)*n_v.v[2];
//      s_v.v[0] = 0.;
//      p_v.v[0] = (4./3.)*p_v.v[1] - (1./3.)*p_v.v[2];
//      q_v.v[0] = 0.;
//      phi_v.v[0] = (4./3.)*phi_v.v[1] - (1./3.)*phi_v.v[2];
//      {
//        int i = nx-1;
//        p_v.v[i] = 0.;
//        q_v.v[i] = 0.;
//        phi_v.v[i] = 0.;
//        s_v.v[0] = 0.;
//        n_v.v[0] = 1;
//      }
//       for(int i = 1; i<nx-1;i++){
//           r_Der_nn = (n_v.v[i+1] - n_v.v[i-1])/(2.*dr[i]);
//           r_Der_ss = (s_v.v[i+1] - s_v.v[i-1])/(2.*dr[i]);
//           r_Der_P =(p_v.v[i+1] - p_v.v[i-1])/(2.*dr[i]);
//           r_Der_Q = (q_v.v[i+1] - q_v.v[i-1])/(2.*dr[i]);
//
//           p_np1[i] = p_v.v[i] +  dt*rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//           q_np1[i] = q_v.v[i] +  dt*rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
//           phi_np1[i] = phi_v.v[i] +  dt*rhs_phi(r[i], n_v.v[i], s_v.v[i], p_v.v[i], q_v.v[i]);
//
//         }
//
//       for(int i = 1; i<nx-1;i++){
//         p_v.v[i] = p_np1[i];
//         q_v.v[i] = q_np1[i];
//         phi_v.v[i] = phi_np1[i];
//       }
//       p_v.check_isfinite(grid.t_evolve);
//       q_v.check_isfinite(grid.t_evolve);
//       phi_v.check_isfinite(grid.t_evolve);
    //==========================================================================
    // Solve_metric_fields solve_metric(grid);
  //   vector<double> p_k1(nx,0),p_k2(nx,0),p_k3(nx,0),p_k4(nx,0);
  //   vector<double> q_k1(nx,0),q_k2(nx,0),q_k3(nx,0),q_k4(nx,0);
  //   vector<double> phi_k1(nx,0),phi_k2(nx,0),phi_k3(nx,0),phi_k4(nx,0);
  //   double r_Der_nn=0., r_Der_ss = 0., r_Der_P = 0., r_Der_Q = 0.;
  //   for(int i = 1; i<nx-1;i++){
  //     if(i==1){
  //       r_Der_nn = Dx_ptp1_4th(n_v.v[i+3], n_v.v[i+2], n_v.v[i+1], n_v.v[i], n_v.v[i-1], dr[i]);
  //       r_Der_ss = Dx_ptp1_4th(s_v.v[i+3], s_v.v[i+2], s_v.v[i+1], s_v.v[i], s_v.v[i-1], dr[i]);
  //       r_Der_P = Dx_ptp1_4th(p_v.v[i+3], p_v.v[i+2], p_v.v[i+1], p_v.v[i], p_v.v[i-1], dr[i]);
  //       r_Der_Q = Dx_ptp1_4th(q_v.v[i+3], q_v.v[i+2], q_v.v[i+1], q_v.v[i], q_v.v[i-1], dr[i]);
  //
  //       p_k1[i] = dt*rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //       q_k1[i] = dt*rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //       phi_k1[i] = dt*rhs_phi(r[i], n_v.v[i], s_v.v[i], p_v.v[i], q_v.v[i]);
  //
  //     }
  //     if(i>1 && i<nx-2){
  //       r_Der_nn = Dx_ptc_4th(n_v.v[i+2], n_v.v[i+1], n_v.v[i-1], n_v.v[i-2], dr[i]);
  //       r_Der_ss = Dx_ptc_4th(s_v.v[i+2], s_v.v[i+1], s_v.v[i-1], s_v.v[i-2], dr[i]);
  //       r_Der_P = Dx_ptc_4th(p_v.v[i+2], p_v.v[i+1], p_v.v[i-1], p_v.v[i-2], dr[i]);
  //       r_Der_Q = Dx_ptc_4th(q_v.v[i+2], q_v.v[i+1], q_v.v[i-1], q_v.v[i-2], dr[i]);
  //
  //       p_k1[i] = dt*rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //       q_k1[i] = dt*rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //       phi_k1[i] = dt*rhs_phi(r[i], n_v.v[i], s_v.v[i], p_v.v[i], q_v.v[i]);
  //
  //     }
  //     if(i==nx-2){
  //       r_Der_nn = Dx_ptm1_4th(n_v.v[i+1], n_v.v[i], n_v.v[i-1], n_v.v[i-2], n_v.v[i-3], dr[i]);
  //       r_Der_ss = Dx_ptm1_4th(s_v.v[i+1], s_v.v[i], s_v.v[i-1], s_v.v[i-2], s_v.v[i-3], dr[i]);
  //       r_Der_P =Dx_ptm1_4th(p_v.v[i+1], p_v.v[i], p_v.v[i-1], p_v.v[i-2], p_v.v[i-3], dr[i]);
  //       r_Der_Q = Dx_ptm1_4th(q_v.v[i+1], q_v.v[i], q_v.v[i-1], q_v.v[i-2], q_v.v[i-3], dr[i]);
  //
  //       p_k1[i] = dt*rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //       q_k1[i] = dt*rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //       phi_k1[i] = dt*rhs_phi(r[i], n_v.v[i], s_v.v[i], p_v.v[i], q_v.v[i]);
  //
  //     }
  //
  //   }
  //   {
  //     Field n_1("n_k1", "even", grid);
  //     Field s_1("s_k1", "odd", grid);
  //     Field p_1("p_k1", "even", grid);
  //     Field q_1("q_k1", "odd", grid);
  //
  //     for(int i = 1; i<nx-1; i++){
  //       n_1.v[i] = n_v.v[i];
  //       s_1.v[i] = s_v.v[i];
  //       p_1.v[i] = p_v.v[i] + 0.5*p_k1[i];
  //       q_1.v[i] = q_v.v[i] + 0.5*q_k1[i];
  //
  //     }
  //     {
  //       int i = nx-1;
  //       n_1.v[i] = n_v.v[i];
  //       s_1.v[i] = s_v.v[i];
  //       p_1.v[i] = 0.;
  //       q_1.v[i] = 0.;
  //     }
  //     {
  //       int i = 0;
  //       n_1.v[i] = make_Dx_zero(n_1.v[i+4], n_1.v[i+3], n_1.v[i+2], n_1.v[i+1]);
  //       s_1.v[i] = 1e-20;
  //       p_1.v[i] = make_Dx_zero(p_1.v[i+4], p_1.v[i+3], p_1.v[i+2], p_1.v[i+1]);
  //       q_1.v[i] = 1e-20;
  //
  //     }
  //     // solve_metric.solve(n_1, s_1, p_1, q_1);
  //
  //     for(int i = 1; i<nx-1;i++){
  //       if(i==1){
  //         r_Der_nn = Dx_ptp1_4th(n_1.v[i+3], n_1.v[i+2], n_1.v[i+1], n_1.v[i], n_1.v[i-1], dr[i]);
  //         r_Der_ss = Dx_ptp1_4th(s_1.v[i+3], s_1.v[i+2], s_1.v[i+1], s_1.v[i], s_1.v[i-1], dr[i]);
  //         r_Der_P = Dx_ptp1_4th(p_1.v[i+3], p_1.v[i+2], p_1.v[i+1], p_1.v[i], p_1.v[i-1], dr[i]);
  //         r_Der_Q = Dx_ptp1_4th(q_1.v[i+3], q_1.v[i+2], q_1.v[i+1], q_1.v[i], q_1.v[i-1], dr[i]);
  //
  //         p_k2[i] = dt*rhs_p(r[i], n_1.v[i], r_Der_nn, s_1.v[i], r_Der_ss, p_1.v[i], r_Der_P, q_1.v[i], r_Der_Q);
  //         q_k2[i] = dt*rhs_q(r[i], n_1.v[i], r_Der_nn, s_1.v[i], r_Der_ss, p_1.v[i], r_Der_P, q_1.v[i], r_Der_Q);
  //         phi_k2[i] = dt*rhs_phi(r[i], n_1.v[i], s_1.v[i], p_1.v[i], q_1.v[i]);
  //
  //       }
  //       if(i>1 && i<nx-2){
  //         r_Der_nn = Dx_ptc_4th(n_1.v[i+2], n_1.v[i+1], n_1.v[i-1], n_1.v[i-2], dr[i]);
  //         r_Der_ss = Dx_ptc_4th(s_1.v[i+2], s_1.v[i+1], s_1.v[i-1], s_1.v[i-2], dr[i]);
  //         r_Der_P = Dx_ptc_4th(p_1.v[i+2], p_1.v[i+1], p_1.v[i-1], p_1.v[i-2], dr[i]);
  //         r_Der_Q = Dx_ptc_4th(q_1.v[i+2], q_1.v[i+1], q_1.v[i-1], q_1.v[i-2], dr[i]);
  //
  //         p_k2[i] = dt*rhs_p(r[i], n_1.v[i], r_Der_nn, s_1.v[i], r_Der_ss, p_1.v[i], r_Der_P, q_1.v[i], r_Der_Q);
  //         q_k2[i] = dt*rhs_q(r[i], n_1.v[i], r_Der_nn, s_1.v[i], r_Der_ss, p_1.v[i], r_Der_P, q_1.v[i], r_Der_Q);
  //         phi_k2[i] = dt*rhs_phi(r[i], n_1.v[i], s_1.v[i], p_1.v[i], q_1.v[i]);
  //
  //       }
  //       if(i==nx-2){
  //         r_Der_nn = Dx_ptm1_4th(n_1.v[i+1], n_1.v[i], n_1.v[i-1], n_1.v[i-2], n_1.v[i-3], dr[i]);
  //         r_Der_ss = Dx_ptm1_4th(s_1.v[i+1], s_1.v[i], s_1.v[i-1], s_1.v[i-2], s_1.v[i-3], dr[i]);
  //         r_Der_P =Dx_ptm1_4th(p_1.v[i+1], p_1.v[i], p_1.v[i-1], p_1.v[i-2], p_1.v[i-3], dr[i]);
  //         r_Der_Q = Dx_ptm1_4th(q_1.v[i+1], q_1.v[i], q_1.v[i-1], q_1.v[i-2], q_1.v[i-3], dr[i]);
  //
  //         p_k2[i] = dt*rhs_p(r[i], n_1.v[i], r_Der_nn, s_1.v[i], r_Der_ss, p_1.v[i], r_Der_P, q_1.v[i], r_Der_Q);
  //         q_k2[i] = dt*rhs_q(r[i], n_1.v[i], r_Der_nn, s_1.v[i], r_Der_ss, p_1.v[i], r_Der_P, q_1.v[i], r_Der_Q);
  //         phi_k2[i] = dt*rhs_phi(r[i], n_1.v[i], s_1.v[i], p_1.v[i], q_1.v[i]);
  //
  //       }
  //
  //     }
  //   }
  //
  //   {
  //     Field n_2("n_k2", "even", grid);
  //     Field s_2("s_k2", "odd", grid);
  //     Field p_2("p_k2", "even", grid);
  //     Field q_2("q_k2", "odd", grid);
  //
  //     for(int i = 1; i<nx-1; i++){
  //       n_2.v[i] = n_v.v[i];
  //       s_2.v[i] = s_v.v[i];
  //       p_2.v[i] = p_v.v[i] + 0.5*p_k2[i];
  //       q_2.v[i] = q_v.v[i] + 0.5*q_k2[i];
  //
  //     }
  //     {
  //       int i = nx-1;
  //       n_2.v[i] = n_v.v[i];
  //       s_2.v[i] = s_v.v[i];
  //       p_2.v[i] = 0.;
  //       q_2.v[i] = 0.;
  //     }
  //     {
  //       int i = 0;
  //       n_2.v[i] = make_Dx_zero(n_2.v[i+4], n_2.v[i+3], n_2.v[i+2], n_2.v[i+1]);
  //       s_2.v[i] = 1e-20;
  //       p_2.v[i] = make_Dx_zero(p_2.v[i+4], p_2.v[i+3], p_2.v[i+2], p_2.v[i+1]);
  //       q_2.v[i] = 1e-20;
  //
  //     }
  //     // solve_metric.solve(n_2, s_2, p_2, q_2);
  //
  //     for(int i = 1; i<nx-1;i++){
  //       if(i==1){
  //         r_Der_nn = Dx_ptp1_4th(n_2.v[i+3], n_2.v[i+2], n_2.v[i+1], n_2.v[i], n_2.v[i-1], dr[i]);
  //         r_Der_ss = Dx_ptp1_4th(s_2.v[i+3], s_2.v[i+2], s_2.v[i+1], s_2.v[i], s_2.v[i-1], dr[i]);
  //         r_Der_P = Dx_ptp1_4th(p_2.v[i+3], p_2.v[i+2], p_2.v[i+1], p_2.v[i], p_2.v[i-1], dr[i]);
  //         r_Der_Q = Dx_ptp1_4th(q_2.v[i+3], q_2.v[i+2], q_2.v[i+1], q_2.v[i], q_2.v[i-1], dr[i]);
  //
  //         p_k3[i] = dt*rhs_p(r[i], n_2.v[i], r_Der_nn, s_2.v[i], r_Der_ss, p_2.v[i], r_Der_P, q_2.v[i], r_Der_Q);
  //         q_k3[i] = dt*rhs_q(r[i], n_2.v[i], r_Der_nn, s_2.v[i], r_Der_ss, p_2.v[i], r_Der_P, q_2.v[i], r_Der_Q);
  //         phi_k3[i] = dt*rhs_phi(r[i], n_2.v[i], s_2.v[i], p_2.v[i], q_2.v[i]);
  //
  //       }
  //       if(i>1 && i<nx-2){
  //         r_Der_nn = Dx_ptc_4th(n_2.v[i+2], n_2.v[i+1], n_2.v[i-1], n_2.v[i-2], dr[i]);
  //         r_Der_ss = Dx_ptc_4th(s_2.v[i+2], s_2.v[i+1], s_2.v[i-1], s_2.v[i-2], dr[i]);
  //         r_Der_P = Dx_ptc_4th(p_2.v[i+2], p_2.v[i+1], p_2.v[i-1], p_2.v[i-2], dr[i]);
  //         r_Der_Q = Dx_ptc_4th(q_2.v[i+2], q_2.v[i+1], q_2.v[i-1], q_2.v[i-2], dr[i]);
  //
  //         p_k3[i] = dt*rhs_p(r[i], n_2.v[i], r_Der_nn, s_2.v[i], r_Der_ss, p_2.v[i], r_Der_P, q_2.v[i], r_Der_Q);
  //         q_k3[i] = dt*rhs_q(r[i], n_2.v[i], r_Der_nn, s_2.v[i], r_Der_ss, p_2.v[i], r_Der_P, q_2.v[i], r_Der_Q);
  //         phi_k3[i] = dt*rhs_phi(r[i], n_2.v[i], s_2.v[i], p_2.v[i], q_2.v[i]);
  //
  //       }
  //       if(i==nx-2){
  //         r_Der_nn = Dx_ptm1_4th(n_2.v[i+1], n_2.v[i], n_2.v[i-1], n_2.v[i-2], n_2.v[i-3], dr[i]);
  //         r_Der_ss = Dx_ptm1_4th(s_2.v[i+1], s_2.v[i], s_2.v[i-1], s_2.v[i-2], s_2.v[i-3], dr[i]);
  //         r_Der_P =Dx_ptm1_4th(p_2.v[i+1], p_2.v[i], p_2.v[i-1], p_2.v[i-2], p_2.v[i-3], dr[i]);
  //         r_Der_Q = Dx_ptm1_4th(q_2.v[i+1], q_2.v[i], q_2.v[i-1], q_2.v[i-2], q_2.v[i-3], dr[i]);
  //
  //         p_k3[i] = dt*rhs_p(r[i], n_2.v[i], r_Der_nn, s_2.v[i], r_Der_ss, p_2.v[i], r_Der_P, q_2.v[i], r_Der_Q);
  //         q_k3[i] = dt*rhs_q(r[i], n_2.v[i], r_Der_nn, s_2.v[i], r_Der_ss, p_2.v[i], r_Der_P, q_2.v[i], r_Der_Q);
  //         phi_k3[i] = dt*rhs_phi(r[i], n_2.v[i], s_2.v[i], p_2.v[i], q_2.v[i]);
  //
  //       }
  //
  //     }
  //
  //
  //
  //
  //   }
  //   {
  //     Field n_3("n_k3", "even", grid);
  //     Field s_3("s_k3", "odd", grid);
  //     Field p_3("p_k3", "even", grid);
  //     Field q_3("q_k3", "odd", grid);
  //
  //     for(int i = 1; i<nx-1; i++){
  //       n_3.v[i] = n_v.v[i];
  //       s_3.v[i] = s_v.v[i];
  //       p_3.v[i] = p_v.v[i] + p_k3[i];
  //       q_3.v[i] = q_v.v[i] + q_k3[i];
  //
  //     }
  //     {
  //       int i = nx-1;
  //       n_3.v[i] = n_v.v[i];
  //       s_3.v[i] = s_v.v[i];
  //       p_3.v[i] = 0.;
  //       q_3.v[i] = 0.;
  //     }
  //     {
  //       int i = 0;
  //       n_3.v[i] = make_Dx_zero(n_3.v[i+4], n_3.v[i+3], n_3.v[i+2], n_3.v[i+1]);
  //       s_3.v[i] = 1e-20;
  //       p_3.v[i] = make_Dx_zero(p_3.v[i+4], p_3.v[i+3], p_3.v[i+2], p_3.v[i+1]);
  //       q_3.v[i] = 1e-20;
  //
  //     }
  //     // solve_metric.solve(n_3, s_3, p_3, q_3);
  //
  //     for(int i = 1; i<nx-1;i++){
  //       if(i==1){
  //         r_Der_nn = Dx_ptp1_4th(n_3.v[i+3], n_3.v[i+2], n_3.v[i+1], n_3.v[i], n_3.v[i-1], dr[i]);
  //         r_Der_ss = Dx_ptp1_4th(s_3.v[i+3], s_3.v[i+2], s_3.v[i+1], s_3.v[i], s_3.v[i-1], dr[i]);
  //         r_Der_P = Dx_ptp1_4th(p_3.v[i+3], p_3.v[i+2], p_3.v[i+1], p_3.v[i], p_3.v[i-1], dr[i]);
  //         r_Der_Q = Dx_ptp1_4th(q_3.v[i+3], q_3.v[i+2], q_3.v[i+1], q_3.v[i], q_3.v[i-1], dr[i]);
  //
  //         p_k4[i] = dt*rhs_p(r[i], n_3.v[i], r_Der_nn, s_3.v[i], r_Der_ss, p_3.v[i], r_Der_P, q_3.v[i], r_Der_Q);
  //         q_k4[i] = dt*rhs_q(r[i], n_3.v[i], r_Der_nn, s_3.v[i], r_Der_ss, p_3.v[i], r_Der_P, q_3.v[i], r_Der_Q);
  //         phi_k4[i] = dt*rhs_phi(r[i], n_3.v[i], s_3.v[i], p_3.v[i], q_3.v[i]);
  //
  //       }
  //       if(i>1 && i<nx-2){
  //         r_Der_nn = Dx_ptc_4th(n_3.v[i+2], n_3.v[i+1], n_3.v[i-1], n_3.v[i-2], dr[i]);
  //         r_Der_ss = Dx_ptc_4th(s_3.v[i+2], s_3.v[i+1], s_3.v[i-1], s_3.v[i-2], dr[i]);
  //         r_Der_P = Dx_ptc_4th(p_3.v[i+2], p_3.v[i+1], p_3.v[i-1], p_3.v[i-2], dr[i]);
  //         r_Der_Q = Dx_ptc_4th(q_3.v[i+2], q_3.v[i+1], q_3.v[i-1], q_3.v[i-2], dr[i]);
  //
  //         p_k4[i] = dt*rhs_p(r[i], n_3.v[i], r_Der_nn, s_3.v[i], r_Der_ss, p_3.v[i], r_Der_P, q_3.v[i], r_Der_Q);
  //         q_k4[i] = dt*rhs_q(r[i], n_3.v[i], r_Der_nn, s_3.v[i], r_Der_ss, p_3.v[i], r_Der_P, q_3.v[i], r_Der_Q);
  //         phi_k4[i] = dt*rhs_phi(r[i], n_3.v[i], s_3.v[i], p_3.v[i], q_3.v[i]);
  //
  //       }
  //       if(i==nx-2){
  //         r_Der_nn = Dx_ptm1_4th(n_3.v[i+1], n_3.v[i], n_3.v[i-1], n_3.v[i-2], n_3.v[i-3], dr[i]);
  //         r_Der_ss = Dx_ptm1_4th(s_3.v[i+1], s_3.v[i], s_3.v[i-1], s_3.v[i-2], s_3.v[i-3], dr[i]);
  //         r_Der_P =Dx_ptm1_4th(p_3.v[i+1], p_3.v[i], p_3.v[i-1], p_3.v[i-2], p_3.v[i-3], dr[i]);
  //         r_Der_Q = Dx_ptm1_4th(q_3.v[i+1], q_3.v[i], q_3.v[i-1], q_3.v[i-2], q_3.v[i-3], dr[i]);
  //
  //         p_k4[i] = dt*rhs_p(r[i], n_3.v[i], r_Der_nn, s_3.v[i], r_Der_ss, p_3.v[i], r_Der_P, q_3.v[i], r_Der_Q);
  //         q_k4[i] = dt*rhs_q(r[i], n_3.v[i], r_Der_nn, s_3.v[i], r_Der_ss, p_3.v[i], r_Der_P, q_3.v[i], r_Der_Q);
  //         phi_k4[i] = dt*rhs_phi(r[i], n_3.v[i], s_3.v[i], p_3.v[i], q_3.v[i]);
  //
  //       }
  //
  //     }
  //
  //   }
  //   for(int i = 1; i< nx-1; i++){
  //     p_v.v[i] = p_v.v[i] + (1./6.)*p_k1[i] + (1./3.)*p_k2[i] + (1./3.)*p_k3[i] + (1./6.)*p_k4[i];
  //     q_v.v[i] = q_v.v[i] + (1./6.)*q_k1[i] + (1./3.)*q_k2[i] + (1./3.)*q_k3[i] + (1./6.)*q_k4[i];
  //     phi_v.v[i] = phi_v.v[i] + (1./6.)*phi_k1[i] + (1./3.)*phi_k2[i] + (1./3.)*phi_k3[i] + (1./6.)*phi_k4[i];
  //
  //   }
  //

  //
  // }


// }
// }
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
//==============================================================================
