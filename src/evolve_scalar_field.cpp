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
  double Q, double r_Der_Q)
{
  return nn*(P + ss*Q);
}
//==============================================================================
double Evolve_scalar_field::rhs_q(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q)
{

  return r_Der_P*nn + r_Der_ss*nn*Q + r_Der_Q*nn*ss + r_Der_nn*(P + Q*ss);
  // return r_Der_P*nn;
}
//==============================================================================
double Evolve_scalar_field::rhs_p(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp)
{

    double Qr = Q/r;
    double ssr = ss/r;
    if((fabs(ss)<1e-6)&& r<10.){
      return Qr*r*r_Der_nn + 2*Qr*nn + r_Der_Q*nn;
    }
    else{
    double num = 0.;
    double denom = 0.;
    {
    double V = 0.;
    double sd0 = 0.;
    double sd1 = 0;
    double sd2 = 0.;
    double sd3 = 0.;
    double sd4 = 0;
    double sd5 = 0;

    V = -1 + 8*Bep*Qr + 8*Bep*P*ssr;
    sd0 = pow(-1 + 8*Bep*Qr,3);
    sd1 = 28*Bep*P*pow(-1 + 8*Bep*Qr,2)*ssr;
    sd2 = 32*pow(Bep,2)*(-1 + 8*Bep*Qr)*(9*pow(P,2) + pow(r,2)*pow(Qr,2))*pow(ssr,2);
    sd3 = 192*pow(Bep,3)*P*(5*pow(P,2) + pow(r,2)*pow(Qr,2))*pow(ssr,3);
    sd4 = 768*pow(Bep,3)*r_Der_Q*(-1 + 8*Bep*Qr)*pow(ssr,4) + 96*pow(Bep,2)*(-1 + 8*Bep*Qr)*(-1 + 8*Bepp*pow(r,2)*pow(Qr,2))*pow(ssr,4);
    sd5 = 6144*pow(Bep,4)*r_Der_Q*P*pow(ssr,5) + 768*pow(Bep,3)*P*(-1 + 8*Bepp*pow(r,2)*pow(Qr,2))*pow(ssr,5);

    denom = V*(sd0 + sd1 + sd2 + sd3 + sd4 + sd5);

    }
    {
      double snm1 = 0.;
      double sn0 = 0.;
      double sn1 = 0.;
      double sn2 = 0.;
      double sn3 = 0.;
      double sn4 = 0.;
      double sn5 =0.;
      double sn6 = 0.;
      double sn7 = 0.;

      snm1 = -0.25*(nn*P*pow(-1 + 8*Bep*Qr,3)*(pow(P,2) - pow(r,2)*pow(Qr,2)))/ssr;
      sn0 = r_Der_Q*nn*pow(-1 + 8*Bep*Qr,2)*(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2)) - nn*pow(-1 + 8*Bep*Qr,2)*(7*Bep*pow(P,4) - 2*Qr + 32*Bep*pow(Qr,2) - 9*Bep*pow(r,2)*pow(P,2)*pow(Qr,2) - 128*pow(Bep,2)*pow(Qr,3) + 2*Bep*pow(r,4)*pow(Qr,4));
      sn1 = r*r_Der_P*nn*pow(-1 + 8*Bep*Qr,3)*(-1 + 12*Bep*Qr)*ssr + 32*r_Der_Q*nn*(-1 + 8*Bep*Qr)*(Bep*P - 16*pow(Bep,2)*P*Qr + 64*pow(Bep,3)*P*pow(Qr,2))*ssr - (nn*(-1 + 8*Bep*Qr)*(3*P + 120*pow(Bep,2)*pow(P,5) - 224*Bep*P*Qr + 3008*pow(Bep,2)*P*pow(Qr,2) - 208*pow(Bep,2)*pow(r,2)*pow(P,3)*pow(Qr,2) - 11264*pow(Bep,3)*P*pow(Qr,3) + 88*pow(Bep,2)*pow(r,4)*P*pow(Qr,4))*ssr)/2.;
      sn2 = 12*Bep*r*r_Der_P*nn*P*pow(-1 + 8*Bep*Qr,2)*(-3 + 44*Bep*Qr)*pow(ssr,2) + 16*Bep*r_Der_Q*nn*(22*Bep*pow(P,2) - 352*pow(Bep,2)*pow(P,2)*Qr - 3*Bep*pow(r,2)*pow(Qr,2) + 1408*pow(Bep,3)*pow(P,2)*pow(Qr,2) + 48*pow(Bep,2)*pow(r,2)*pow(Qr,3) - 192*pow(Bep,3)*pow(r,2)*pow(Qr,4))*pow(ssr,2) + 2*Bep*nn*(-32*pow(P,2) - 16*Bepp*pow(P,4) - 80*pow(Bep,2)*pow(P,6) + 1296*Bep*pow(P,2)*Qr + 256*Bep*Bepp*pow(P,4)*Qr + 3*pow(r,2)*pow(Qr,2) - 14592*pow(Bep,2)*pow(P,2)*pow(Qr,2) + 40*Bepp*pow(r,2)*pow(P,2)*pow(Qr,2) - 1024*pow(Bep,2)*Bepp*pow(P,4)*pow(Qr,2) + 184*pow(Bep,2)*pow(r,2)*pow(P,4)*pow(Qr,2) - 64*Bep*pow(r,2)*pow(Qr,3) + 50176*pow(Bep,3)*pow(P,2)*pow(Qr,3) - 640*Bep*Bepp*pow(r,2)*pow(P,2)*pow(Qr,3) + 448*pow(Bep,2)*pow(r,2)*pow(Qr,4) - 24*Bepp*pow(r,4)*pow(Qr,4) + 2560*pow(Bep,2)*Bepp*pow(r,2)*pow(P,2)*pow(Qr,4) - 112*pow(Bep,2)*pow(r,4)*pow(P,2)*pow(Qr,4) - 1024*pow(Bep,3)*pow(r,2)*pow(Qr,5) + 384*Bep*Bepp*pow(r,4)*pow(Qr,5) - 1536*pow(Bep,2)*Bepp*pow(r,4)*pow(Qr,6) + 8*pow(Bep,2)*pow(r,6)*pow(Qr,6))*pow(ssr,2);

      sn3 = 32*pow(Bep,2)*r_Der_P*nn*(-1 + 8*Bep*Qr)*(-16*r*pow(P,2) + 232*Bep*r*pow(P,2)*Qr - pow(r,3)*pow(Qr,2))*pow(ssr,3) - 128*pow(Bep,2)*r_Der_Q*nn*(-1 + 8*Bep*Qr)*(-12*Bep*pow(P,3) + 7*Bep*pow(r,2)*P*pow(Qr,2))*pow(ssr,3) - 8*pow(Bep,2)*nn*(-1 + 8*Bep*Qr)*(121*pow(P,3) + 56*Bepp*pow(P,5) - 1776*Bep*pow(P,3)*Qr - 21*pow(r,2)*P*pow(Qr,2) - 200*Bepp*pow(r,2)*pow(P,3)*pow(Qr,2) + 112*Bep*pow(r,2)*P*pow(Qr,3) + 144*Bepp*pow(r,4)*P*pow(Qr,4))*pow(ssr,3);

      sn4 = 768*Bep*pow(r_Der_P,2)*nn*(pow(Bep,2) - 16*pow(Bep,3)*Qr + 64*pow(Bep,4)*pow(Qr,2))*pow(ssr,4) + 64*Bep*r_Der_P*nn*(-51*pow(Bep,2)*r*pow(P,3) + 24*Bep*Bepp*r*P*Qr + 612*pow(Bep,3)*r*pow(P,3)*Qr - 384*pow(Bep,2)*Bepp*r*P*pow(Qr,2) - 7*pow(Bep,2)*pow(r,3)*P*pow(Qr,2) + 1536*pow(Bep,3)*Bepp*r*P*pow(Qr,3) + 20*pow(Bep,3)*pow(r,3)*P*pow(Qr,3))*pow(ssr,4) + 32*Bep*r_Der_Q*nn*(-3*Bep - 24*Bep*Bepp*pow(P,2) + 64*pow(Bep,3)*pow(P,4) + 48*pow(Bep,2)*Qr + 384*pow(Bep,2)*Bepp*pow(P,2)*Qr - 192*pow(Bep,3)*pow(Qr,2) - 1536*pow(Bep,3)*Bepp*pow(P,2)*pow(Qr,2) - 120*pow(Bep,3)*pow(r,2)*pow(P,2)*pow(Qr,2) + 8*pow(Bep,3)*pow(r,4)*pow(Qr,4))*pow(ssr,4) - 4*Bep*nn*(-3 - 24*Bepp*pow(P,2) + 1584*pow(Bep,2)*pow(P,4) + 384*pow(Bep,2)*Bepp*pow(P,6) + 48*Bep*Qr + 384*Bep*Bepp*pow(P,2)*Qr - 16384*pow(Bep,3)*pow(P,4)*Qr - 192*pow(Bep,2)*pow(Qr,2) + 24*Bepp*pow(r,2)*pow(Qr,2) - 1536*pow(Bep,2)*Bepp*pow(P,2)*pow(Qr,2) - 328*pow(Bep,2)*pow(r,2)*pow(P,2)*pow(Qr,2) - 1856*pow(Bep,2)*Bepp*pow(r,2)*pow(P,4)*pow(Qr,2) - 384*Bep*Bepp*pow(r,2)*pow(Qr,3) + 2048*pow(Bep,3)*pow(r,2)*pow(P,2)*pow(Qr,3) + 1536*pow(Bep,2)*Bepp*pow(r,2)*pow(Qr,4) + 8*pow(Bep,2)*pow(r,4)*pow(Qr,4) + 1536*pow(Bep,2)*Bepp*pow(r,4)*pow(P,2)*pow(Qr,4) - 64*pow(Bep,2)*Bepp*pow(r,6)*pow(Qr,6))*pow(ssr,4);

      sn5 = 15360*pow(Bep,4)*pow(r_Der_P,2)*nn*P*(-1 + 8*Bep*Qr)*pow(ssr,5) - 1536*pow(Bep,2)*r_Der_Q*nn*(-(Bep*P) - 8*Bep*Bepp*pow(P,3) + 8*pow(Bep,2)*P*Qr + 64*pow(Bep,2)*Bepp*pow(P,3)*Qr - 2*Bep*Bepp*pow(r,2)*P*pow(Qr,2) + 16*pow(Bep,2)*Bepp*pow(r,2)*P*pow(Qr,3))*pow(ssr,5) + 192*pow(Bep,2)*nn*(-P - 8*Bepp*pow(P,3) + 80*pow(Bep,2)*pow(P,5) + 8*Bep*P*Qr + 64*Bep*Bepp*pow(P,3)*Qr + 6*Bepp*pow(r,2)*P*pow(Qr,2) - 16*pow(Bep,2)*pow(r,2)*pow(P,3)*pow(Qr,2) - 16*pow(Bepp,2)*pow(r,2)*pow(P,3)*pow(Qr,2) - 48*Bep*Bepp*pow(r,2)*P*pow(Qr,3) + 128*Bep*pow(Bepp,2)*pow(r,2)*pow(P,3)*pow(Qr,3) + 16*pow(Bepp,2)*pow(r,4)*P*pow(Qr,4) - 128*Bep*pow(Bepp,2)*pow(r,4)*P*pow(Qr,5))*pow(ssr,5) + r_Der_P*(768*pow(Bep,3)*r*r_Der_Q*nn*(1 - 12*Bep*Qr + 32*pow(Bep,2)*pow(Qr,2))*pow(ssr,5) + 96*pow(Bep,2)*nn*(-r + 80*pow(Bep,2)*r*pow(P,4) + 12*Bep*r*Qr - 320*Bep*Bepp*r*pow(P,2)*Qr - 32*pow(Bep,2)*r*pow(Qr,2) + 8*Bepp*pow(r,3)*pow(Qr,2) + 2560*pow(Bep,2)*Bepp*r*pow(P,2)*pow(Qr,2) + 16*pow(Bep,2)*pow(r,3)*pow(P,2)*pow(Qr,2) - 96*Bep*Bepp*pow(r,3)*pow(Qr,3) + 256*pow(Bep,2)*Bepp*pow(r,3)*pow(Qr,4))*pow(ssr,5));

      sn6 = 73728*pow(Bep,5)*pow(r_Der_P,2)*nn*pow(P,2)*pow(ssr,6) - 6144*pow(Bep,3)*r_Der_Q*nn*P*(Bep*P + 8*Bep*Bepp*pow(P,3) + 4*Bep*Bepp*pow(r,2)*P*pow(Qr,2))*pow(ssr,6) + 768*pow(Bep,3)*nn*P*(P + 8*Bepp*pow(P,3) - 4*Bepp*pow(r,2)*P*pow(Qr,2) + 32*pow(Bepp,2)*pow(r,2)*pow(P,3)*pow(Qr,2) - 32*pow(Bepp,2)*pow(r,4)*P*pow(Qr,4))*pow(ssr,6) + r_Der_P*(12288*pow(Bep,4)*r*r_Der_Q*nn*P*(-1 + 6*Bep*Qr)*pow(ssr,6) + 1536*pow(Bep,3)*nn*P*(r - 6*Bep*r*Qr + 96*Bep*Bepp*r*pow(P,2)*Qr - 8*Bepp*pow(r,3)*pow(Qr,2) + 48*Bep*Bepp*pow(r,3)*pow(Qr,3))*pow(ssr,6));

      sn7 = r_Der_P*(49152*pow(Bep,5)*r*r_Der_Q*nn*pow(P,2)*pow(ssr,7) + 6144*pow(Bep,4)*r*nn*pow(P,2)*(-1 + 8*Bepp*pow(r,2)*pow(Qr,2))*pow(ssr,7));

      num = snm1 + sn0 + sn1 + sn2 + sn3 + sn4 + sn5 + sn6 + sn7 ;
    }
    double ans = num;
    ans/= denom;

    return ans;
  }
}
//==============================================================================
double Evolve_scalar_field::rhs_s_free(double r,
  double nn, double r_Der_nn,
  double ss, double r_Der_ss,
  double P, double r_Der_P,
  double Q, double r_Der_Q,
  double Bep, double Bepp)
{

  double Qr = Q/r;
  double ssr = ss/r;
  double ans = 0.;

  // return
  // ((pow(Qr,2)*pow(r,3)*nn)/4. - 2*Bep*pow(Qr,3)*pow(r,3)*nn + (r*pow(ssr,2)*nn)/2. + 4*Bep*Qr*r*pow(ssr,2)*nn - 64*pow(Bep,2)*pow(Qr,2)*r*pow(ssr,2)*nn + 16*pow(Bep,2)*pow(Qr,2)*pow(r,3)*pow(ssr,4)*nn + (Qr*r*nn*P)/(2.*ssr) - (4*Bep*pow(Qr,2)*r*nn*P)/ssr + 4*Bepp*Qr*r*ssr*nn*P - 32*Bep*Bepp*pow(Qr,2)*r*ssr*nn*P - 2*Bep*pow(Qr,2)*pow(r,3)*ssr*nn*P + 4*Bep*r*pow(ssr,3)*nn*P - 128*pow(Bep,2)*Qr*r*pow(ssr,3)*nn*P + (r*nn*pow(P,2))/4. - 6*Bep*Qr*r*nn*pow(P,2) + 4*Bepp*r*pow(ssr,2)*nn*pow(P,2) - 64*Bep*Bepp*Qr*r*pow(ssr,2)*nn*pow(P,2) - 80*pow(Bep,2)*r*pow(ssr,4)*nn*pow(P,2) - 2*Bep*r*ssr*nn*pow(P,3) - 32*Bep*Bepp*r*pow(ssr,3)*nn*pow(P,3) - pow(r_Der_nn,2)*((256*pow(Bep,3)*Qr*pow(r,3)*pow(ssr,6))/nn + (256*pow(Bep,3)*r*pow(ssr,5)*P)/nn) - r_Der_P*(-4*Bep*ssr*nn + 32*pow(Bep,2)*Qr*ssr*nn + 32*pow(Bep,2)*pow(ssr,2)*nn*P) - r_Der_Q*(-4*Bep*r*pow(ssr,2)*nn + 32*pow(Bep,2)*Qr*r*pow(ssr,2)*nn + 32*pow(Bep,2)*r*pow(ssr,3)*nn*P) - pow(r_Der_ss,2)*(-128*pow(Bep,2)*r*pow(ssr,4)*nn + 768*pow(Bep,3)*Qr*r*pow(ssr,4)*nn + 768*pow(Bep,3)*r*pow(ssr,5)*nn*P) - r_Der_ss*(-(r*ssr*nn) + 12*Bep*Qr*r*ssr*nn - 32*pow(Bep,2)*pow(Qr,2)*r*ssr*nn - 16*pow(Bep,2)*pow(Qr,2)*pow(r,3)*pow(ssr,3)*nn - 256*pow(Bep,3)*r_Der_P*pow(ssr,4)*nn + 32*pow(Bep,2)*r*pow(ssr,5)*nn - 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,3)*pow(ssr,5)*nn - 256*pow(Bep,3)*r*r_Der_Q*pow(ssr,5)*nn + 12*Bep*r*pow(ssr,2)*nn*P - 96*pow(Bep,2)*Qr*r*pow(ssr,2)*nn*P - 256*pow(Bep,2)*Bepp*Qr*r*pow(ssr,4)*nn*P - 48*pow(Bep,2)*r*pow(ssr,3)*nn*pow(P,2)) - r_Der_nn*(-4*Bep*Qr*pow(r,2)*pow(ssr,2) + 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,2) - 32*pow(Bep,2)*pow(ssr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4) - 16*pow(Bep,2)*pow(Qr,2)*pow(r,4)*pow(ssr,4) + 256*pow(Bep,3)*r_Der_Q*pow(ssr,4) + 256*pow(Bep,3)*r*r_Der_P*pow(ssr,5) - 4*Bep*pow(r,2)*pow(ssr,3)*P + 32*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,3)*P + 256*pow(Bep,2)*Bepp*Qr*pow(r,2)*pow(ssr,5)*P + 16*pow(Bep,2)*pow(r,2)*pow(ssr,4)*pow(P,2) + r_Der_ss*(64*pow(Bep,2)*pow(ssr,3) - 512*pow(Bep,3)*Qr*pow(ssr,3) - 128*pow(Bep,2)*pow(r,2)*pow(ssr,5) + 1024*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,5) - 256*pow(Bep,3)*pow(ssr,4)*P + 768*pow(Bep,3)*pow(r,2)*pow(ssr,6)*P)))/(1 - 16*Bep*Qr + 64*pow(Bep,2)*pow(Qr,2) - 32*pow(Bep,2)*pow(ssr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4) + 256*pow(Bep,3)*r_Der_Q*pow(ssr,4) - 16*Bep*ssr*P + 128*pow(Bep,2)*Qr*ssr*P + 64*pow(Bep,2)*pow(ssr,2)*pow(P,2) + r_Der_ss*(128*pow(Bep,2)*pow(ssr,3) - 1024*pow(Bep,3)*Qr*pow(ssr,3) - 768*pow(Bep,3)*pow(ssr,4)*P) + r_Der_nn*((128*pow(Bep,2)*r*pow(ssr,4))/nn - (1024*pow(Bep,3)*Qr*r*pow(ssr,4))/nn - (768*pow(Bep,3)*r*pow(ssr,5)*P)/nn))
  // ;
  ans = (-12288*pow(Bep,4)*pow(r,6)*pow(r_Der_P,2)*pow(ssr,7)*pow(nn,2))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) + (12288*pow(r,6)*pow(r_Der_Q,2)*pow(nn,2)*(-(pow(Bep,4)*pow(ssr,7)) + 12*pow(Bep,5)*Qr*pow(ssr,7) - 32*pow(Bep,6)*pow(Qr,2)*pow(ssr,7) - 4*pow(Bep,5)*Qr*pow(r,2)*pow(ssr,9) + 16*pow(Bep,6)*pow(Qr,2)*pow(r,2)*pow(ssr,9) + 16*pow(Bep,5)*pow(ssr,8)*P - 96*pow(Bep,6)*Qr*pow(ssr,8)*P + 32*pow(Bep,6)*Qr*pow(r,2)*pow(ssr,10)*P - 64*pow(Bep,6)*pow(ssr,9)*pow(P,2)))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*pow(-1 + 8*Bep*Qr + 12*Bep*ssr*P,2)) + (32*pow(r,6)*r_Der_Q*pow(nn,2)*(-(Bep*pow(ssr,3)) + 30*pow(Bep,2)*Qr*pow(ssr,3) - 336*pow(Bep,3)*pow(Qr,2)*pow(ssr,3) + 1664*pow(Bep,4)*pow(Qr,3)*pow(ssr,3) - 3072*pow(Bep,5)*pow(Qr,4)*pow(ssr,3) - 2*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,5) + 224*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(ssr,5) - 768*pow(Bep,5)*pow(Qr,4)*pow(r,2)*pow(ssr,5) + 96*pow(Bep,3)*pow(ssr,7) - 1152*pow(Bep,4)*Qr*pow(ssr,7) + 3072*pow(Bep,5)*pow(Qr,2)*pow(ssr,7) - 768*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,7) + 9216*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,7) - 24576*pow(Bep,5)*Bepp*pow(Qr,4)*pow(r,2)*pow(ssr,7) - 160*pow(Bep,4)*pow(Qr,3)*pow(r,4)*pow(ssr,7) + 512*pow(Bep,5)*pow(Qr,4)*pow(r,4)*pow(ssr,7) + 384*pow(Bep,4)*Qr*pow(r,2)*pow(ssr,9) - 1536*pow(Bep,5)*pow(Qr,2)*pow(r,2)*pow(ssr,9) - 3072*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,9) + 12288*pow(Bep,5)*Bepp*pow(Qr,4)*pow(r,4)*pow(ssr,9) + 36*pow(Bep,2)*pow(ssr,4)*P - 856*pow(Bep,3)*Qr*pow(ssr,4)*P + 6784*pow(Bep,4)*pow(Qr,2)*pow(ssr,4)*P - 17920*pow(Bep,5)*pow(Qr,3)*pow(ssr,4)*P - 384*pow(Bep,3)*Bepp*Qr*pow(ssr,6)*P + 6144*pow(Bep,4)*Bepp*pow(Qr,2)*pow(ssr,6)*P - 24576*pow(Bep,5)*Bepp*pow(Qr,3)*pow(ssr,6)*P + 8*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,6)*P + 384*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,6)*P - 2432*pow(Bep,5)*pow(Qr,3)*pow(r,2)*pow(ssr,6)*P - 1536*pow(Bep,4)*pow(ssr,8)*P + 9216*pow(Bep,5)*Qr*pow(ssr,8)*P - 384*pow(Bep,3)*Bepp*Qr*pow(r,2)*pow(ssr,8)*P + 15360*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,8)*P - 73728*pow(Bep,5)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,8)*P + 1152*pow(Bep,5)*pow(Qr,3)*pow(r,4)*pow(ssr,8)*P - 3072*pow(Bep,5)*Qr*pow(r,2)*pow(ssr,10)*P + 24576*pow(Bep,5)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,10)*P - 528*pow(Bep,3)*pow(ssr,5)*pow(P,2) + 8736*pow(Bep,4)*Qr*pow(ssr,5)*pow(P,2) - 36096*pow(Bep,5)*pow(Qr,2)*pow(ssr,5)*pow(P,2) + 7680*pow(Bep,4)*Bepp*Qr*pow(ssr,7)*pow(P,2) - 61440*pow(Bep,5)*Bepp*pow(Qr,2)*pow(ssr,7)*pow(P,2) + 288*pow(Bep,4)*Qr*pow(r,2)*pow(ssr,7)*pow(P,2) - 3072*pow(Bep,5)*pow(Qr,2)*pow(r,2)*pow(ssr,7)*pow(P,2) + 6144*pow(Bep,5)*pow(ssr,9)*pow(P,2) + 7680*pow(Bep,4)*Bepp*Qr*pow(r,2)*pow(ssr,9)*pow(P,2) - 73728*pow(Bep,5)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,9)*pow(P,2) + 3584*pow(Bep,4)*pow(ssr,6)*pow(P,3) - 30592*pow(Bep,5)*Qr*pow(ssr,6)*pow(P,3) - 36864*pow(Bep,5)*Bepp*Qr*pow(ssr,8)*pow(P,3) - 1920*pow(Bep,5)*Qr*pow(r,2)*pow(ssr,8)*pow(P,3) - 36864*pow(Bep,5)*Bepp*Qr*pow(r,2)*pow(ssr,10)*pow(P,3) - 9216*pow(Bep,5)*pow(ssr,7)*pow(P,4)))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*pow(-1 + 8*Bep*Qr + 12*Bep*ssr*P,2)) + (2*pow(r,6)*pow(nn,2)*(-(pow(Qr,2)*pow(r,2)*ssr) + 30*Bep*pow(Qr,3)*pow(r,2)*ssr - 336*pow(Bep,2)*pow(Qr,4)*pow(r,2)*ssr + 1664*pow(Bep,3)*pow(Qr,5)*pow(r,2)*ssr - 3072*pow(Bep,4)*pow(Qr,6)*pow(r,2)*ssr - 12*Bep*Qr*pow(ssr,3) + 416*pow(Bep,2)*pow(Qr,2)*pow(ssr,3) - 5376*pow(Bep,3)*pow(Qr,3)*pow(ssr,3) + 30720*pow(Bep,4)*pow(Qr,4)*pow(ssr,3) - 65536*pow(Bep,5)*pow(Qr,5)*pow(ssr,3) - 8*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,3) + 224*Bep*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,3) - 2304*pow(Bep,2)*Bepp*pow(Qr,4)*pow(r,2)*pow(ssr,3) + 10240*pow(Bep,3)*Bepp*pow(Qr,5)*pow(r,2)*pow(ssr,3) - 16384*pow(Bep,4)*Bepp*pow(Qr,6)*pow(r,2)*pow(ssr,3) - 2*Bep*pow(Qr,3)*pow(r,4)*pow(ssr,3) + 24*pow(Bep,2)*pow(Qr,4)*pow(r,4)*pow(ssr,3) - 64*pow(Bep,3)*pow(Qr,5)*pow(r,4)*pow(ssr,3) + 4*Bep*Qr*pow(r,2)*pow(ssr,5) - 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,5) + 320*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(ssr,5) - 4608*pow(Bep,4)*pow(Qr,4)*pow(r,2)*pow(ssr,5) + 16384*pow(Bep,5)*pow(Qr,5)*pow(r,2)*pow(ssr,5) - 32*Bep*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,5) + 3584*pow(Bep,3)*Bepp*pow(Qr,5)*pow(r,4)*pow(ssr,5) - 12288*pow(Bep,4)*Bepp*pow(Qr,6)*pow(r,4)*pow(ssr,5) - 64*pow(Bep,3)*pow(Qr,5)*pow(r,6)*pow(ssr,5) + 128*pow(Bep,4)*pow(Qr,6)*pow(r,6)*pow(ssr,5) - 96*pow(Bep,2)*pow(ssr,7) + 1152*pow(Bep,3)*Qr*pow(ssr,7) - 3072*pow(Bep,4)*pow(Qr,2)*pow(ssr,7) + 1536*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,7) - 18432*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,7) + 49152*pow(Bep,4)*Bepp*pow(Qr,4)*pow(r,2)*pow(ssr,7) + 320*pow(Bep,3)*pow(Qr,3)*pow(r,4)*pow(ssr,7) - 1024*pow(Bep,4)*pow(Qr,4)*pow(r,4)*pow(ssr,7) - 6144*pow(Bep,2)*pow(Bepp,2)*pow(Qr,4)*pow(r,4)*pow(ssr,7) + 73728*pow(Bep,3)*pow(Bepp,2)*pow(Qr,5)*pow(r,4)*pow(ssr,7) - 196608*pow(Bep,4)*pow(Bepp,2)*pow(Qr,6)*pow(r,4)*pow(ssr,7) - 2560*pow(Bep,3)*Bepp*pow(Qr,5)*pow(r,6)*pow(ssr,7) + 8192*pow(Bep,4)*Bepp*pow(Qr,6)*pow(r,6)*pow(ssr,7) - 384*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,9) + 1536*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,9) + 6144*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,9) - 24576*pow(Bep,4)*Bepp*pow(Qr,4)*pow(r,4)*pow(ssr,9) - 24576*pow(Bep,3)*pow(Bepp,2)*pow(Qr,5)*pow(r,6)*pow(ssr,9) + 98304*pow(Bep,4)*pow(Bepp,2)*pow(Qr,6)*pow(r,6)*pow(ssr,9) - Qr*P + 32*Bep*pow(Qr,2)*P - 384*pow(Bep,2)*pow(Qr,3)*P + 2048*pow(Bep,3)*pow(Qr,4)*P - 4096*pow(Bep,4)*pow(Qr,5)*P - 8*Bepp*Qr*pow(ssr,2)*P + 256*Bep*Bepp*pow(Qr,2)*pow(ssr,2)*P - 3072*pow(Bep,2)*Bepp*pow(Qr,3)*pow(ssr,2)*P + 16384*pow(Bep,3)*Bepp*pow(Qr,4)*pow(ssr,2)*P - 32768*pow(Bep,4)*Bepp*pow(Qr,5)*pow(ssr,2)*P - Qr*pow(r,2)*pow(ssr,2)*P + 68*Bep*pow(Qr,2)*pow(r,2)*pow(ssr,2)*P - 1224*pow(Bep,2)*pow(Qr,3)*pow(r,2)*pow(ssr,2)*P + 8576*pow(Bep,3)*pow(Qr,4)*pow(r,2)*pow(ssr,2)*P - 20992*pow(Bep,4)*pow(Qr,5)*pow(r,2)*pow(ssr,2)*P - 8*Bep*pow(ssr,4)*P + 944*pow(Bep,2)*Qr*pow(ssr,4)*P - 19712*pow(Bep,3)*pow(Qr,2)*pow(ssr,4)*P + 150528*pow(Bep,4)*pow(Qr,3)*pow(ssr,4)*P - 393216*pow(Bep,5)*pow(Qr,4)*pow(ssr,4)*P - 8*Bepp*Qr*pow(r,2)*pow(ssr,4)*P + 512*Bep*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4)*P - 9344*pow(Bep,2)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,4)*P + 67584*pow(Bep,3)*Bepp*pow(Qr,4)*pow(r,2)*pow(ssr,4)*P - 172032*pow(Bep,4)*Bepp*pow(Qr,5)*pow(r,2)*pow(ssr,4)*P + 24*pow(Bep,2)*pow(Qr,3)*pow(r,4)*pow(ssr,4)*P - 256*pow(Bep,3)*pow(Qr,4)*pow(r,4)*pow(ssr,4)*P + 896*pow(Bep,4)*pow(Qr,5)*pow(r,4)*pow(ssr,4)*P + 768*pow(Bep,2)*Bepp*Qr*pow(ssr,6)*P - 12288*pow(Bep,3)*Bepp*pow(Qr,2)*pow(ssr,6)*P + 49152*pow(Bep,4)*Bepp*pow(Qr,3)*pow(ssr,6)*P - 16*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,6)*P + 256*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(ssr,6)*P - 11520*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(ssr,6)*P - 6144*pow(Bep,2)*pow(Bepp,2)*pow(Qr,3)*pow(r,2)*pow(ssr,6)*P + 65536*pow(Bep,5)*pow(Qr,4)*pow(r,2)*pow(ssr,6)*P + 98304*pow(Bep,3)*pow(Bepp,2)*pow(Qr,4)*pow(r,2)*pow(ssr,6)*P - 393216*pow(Bep,4)*pow(Bepp,2)*pow(Qr,5)*pow(r,2)*pow(ssr,6)*P - 128*pow(Bep,2)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,6)*P + 7168*pow(Bep,3)*Bepp*pow(Qr,4)*pow(r,4)*pow(ssr,6)*P - 30720*pow(Bep,4)*Bepp*pow(Qr,5)*pow(r,4)*pow(ssr,6)*P + 384*pow(Bep,4)*pow(Qr,5)*pow(r,6)*pow(ssr,6)*P + 1536*pow(Bep,3)*pow(ssr,8)*P - 9216*pow(Bep,4)*Qr*pow(ssr,8)*P + 768*pow(Bep,2)*Bepp*Qr*pow(r,2)*pow(ssr,8)*P - 30720*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,8)*P + 147456*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,8)*P - 2304*pow(Bep,4)*pow(Qr,3)*pow(r,4)*pow(ssr,8)*P - 6144*pow(Bep,2)*pow(Bepp,2)*pow(Qr,3)*pow(r,4)*pow(ssr,8)*P + 147456*pow(Bep,3)*pow(Bepp,2)*pow(Qr,4)*pow(r,4)*pow(ssr,8)*P - 589824*pow(Bep,4)*pow(Bepp,2)*pow(Qr,5)*pow(r,4)*pow(ssr,8)*P + 18432*pow(Bep,4)*Bepp*pow(Qr,5)*pow(r,6)*pow(ssr,8)*P + 3072*pow(Bep,4)*Qr*pow(r,2)*pow(ssr,10)*P - 49152*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,10)*P + 196608*pow(Bep,4)*pow(Bepp,2)*pow(Qr,5)*pow(r,6)*pow(ssr,10)*P - ssr*pow(P,2) + 70*Bep*Qr*ssr*pow(P,2) - 1296*pow(Bep,2)*pow(Qr,2)*ssr*pow(P,2) + 9344*pow(Bep,3)*pow(Qr,3)*ssr*pow(P,2) - 23552*pow(Bep,4)*pow(Qr,4)*ssr*pow(P,2) - 8*Bepp*pow(ssr,3)*pow(P,2) + 576*Bep*Bepp*Qr*pow(ssr,3)*pow(P,2) - 10752*pow(Bep,2)*Bepp*pow(Qr,2)*pow(ssr,3)*pow(P,2) + 77824*pow(Bep,3)*Bepp*pow(Qr,3)*pow(ssr,3)*pow(P,2) - 196608*pow(Bep,4)*Bepp*pow(Qr,4)*pow(ssr,3)*pow(P,2) + 38*Bep*Qr*pow(r,2)*pow(ssr,3)*pow(P,2) - 1472*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,3)*pow(P,2) + 16192*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(ssr,3)*pow(P,2) - 54784*pow(Bep,4)*pow(Qr,4)*pow(r,2)*pow(ssr,3)*pow(P,2) + 544*pow(Bep,2)*pow(ssr,5)*pow(P,2) - 24128*pow(Bep,3)*Qr*pow(ssr,5)*pow(P,2) + 276992*pow(Bep,4)*pow(Qr,2)*pow(ssr,5)*pow(P,2) - 950272*pow(Bep,5)*pow(Qr,3)*pow(ssr,5)*pow(P,2) + 320*Bep*Bepp*Qr*pow(r,2)*pow(ssr,5)*pow(P,2) - 12416*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,5)*pow(P,2) + 146432*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,5)*pow(P,2) - 540672*pow(Bep,4)*Bepp*pow(Qr,4)*pow(r,2)*pow(ssr,5)*pow(P,2) - 64*pow(Bep,3)*pow(Qr,3)*pow(r,4)*pow(ssr,5)*pow(P,2) + 1280*pow(Bep,4)*pow(Qr,4)*pow(r,4)*pow(ssr,5)*pow(P,2) - 15360*pow(Bep,3)*Bepp*Qr*pow(ssr,7)*pow(P,2) + 122880*pow(Bep,4)*Bepp*pow(Qr,2)*pow(ssr,7)*pow(P,2) - 576*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,7)*pow(P,2) - 4608*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,7)*pow(P,2) - 6144*pow(Bep,2)*pow(Bepp,2)*pow(Qr,2)*pow(r,2)*pow(ssr,7)*pow(P,2) + 86016*pow(Bep,5)*pow(Qr,3)*pow(r,2)*pow(ssr,7)*pow(P,2) + 221184*pow(Bep,3)*pow(Bepp,2)*pow(Qr,3)*pow(r,2)*pow(ssr,7)*pow(P,2) - 1376256*pow(Bep,4)*pow(Bepp,2)*pow(Qr,4)*pow(r,2)*pow(ssr,7)*pow(P,2) + 9216*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,7)*pow(P,2) - 49152*pow(Bep,4)*Bepp*pow(Qr,4)*pow(r,4)*pow(ssr,7)*pow(P,2) - 6144*pow(Bep,4)*pow(ssr,9)*pow(P,2) - 15360*pow(Bep,3)*Bepp*Qr*pow(r,2)*pow(ssr,9)*pow(P,2) + 147456*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,9)*pow(P,2) + 122880*pow(Bep,3)*pow(Bepp,2)*pow(Qr,3)*pow(r,4)*pow(ssr,9)*pow(P,2) - 786432*pow(Bep,4)*pow(Bepp,2)*pow(Qr,4)*pow(r,4)*pow(ssr,9)*pow(P,2) + 36*Bep*pow(ssr,2)*pow(P,3) - 1432*pow(Bep,2)*Qr*pow(ssr,2)*pow(P,3) + 16000*pow(Bep,3)*pow(Qr,2)*pow(ssr,2)*pow(P,3) - 54784*pow(Bep,4)*pow(Qr,3)*pow(ssr,2)*pow(P,3) + 320*Bep*Bepp*pow(ssr,4)*pow(P,3) - 12672*pow(Bep,2)*Bepp*Qr*pow(ssr,4)*pow(P,3) + 141312*pow(Bep,3)*Bepp*pow(Qr,2)*pow(ssr,4)*pow(P,3) - 483328*pow(Bep,4)*Bepp*pow(Qr,3)*pow(ssr,4)*pow(P,3) - 568*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,4)*pow(P,3) + 13184*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(ssr,4)*pow(P,3) - 69120*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(ssr,4)*pow(P,3) - 9984*pow(Bep,3)*pow(ssr,6)*pow(P,3) + 229120*pow(Bep,4)*Qr*pow(ssr,6)*pow(P,3) - 1163264*pow(Bep,5)*pow(Qr,2)*pow(ssr,6)*pow(P,3) - 4992*pow(Bep,2)*Bepp*Qr*pow(r,2)*pow(ssr,6)*pow(P,3) + 128000*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,6)*pow(P,3) - 753664*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,6)*pow(P,3) + 73728*pow(Bep,4)*Bepp*Qr*pow(ssr,8)*pow(P,3) + 3840*pow(Bep,4)*Qr*pow(r,2)*pow(ssr,8)*pow(P,3) + 36864*pow(Bep,5)*pow(Qr,2)*pow(r,2)*pow(ssr,8)*pow(P,3) + 147456*pow(Bep,3)*pow(Bepp,2)*pow(Qr,2)*pow(r,2)*pow(ssr,8)*pow(P,3) - 1769472*pow(Bep,4)*pow(Bepp,2)*pow(Qr,3)*pow(r,2)*pow(ssr,8)*pow(P,3) - 49152*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,4)*pow(ssr,8)*pow(P,3) + 73728*pow(Bep,4)*Bepp*Qr*pow(r,2)*pow(ssr,10)*pow(P,3) - 589824*pow(Bep,4)*pow(Bepp,2)*pow(Qr,3)*pow(r,4)*pow(ssr,10)*pow(P,3) - 504*pow(Bep,2)*pow(ssr,3)*pow(P,4) + 11904*pow(Bep,3)*Qr*pow(ssr,3)*pow(P,4) - 62976*pow(Bep,4)*pow(Qr,2)*pow(ssr,3)*pow(P,4) - 4736*pow(Bep,2)*Bepp*pow(ssr,5)*pow(P,4) + 111104*pow(Bep,3)*Bepp*Qr*pow(ssr,5)*pow(P,4) - 585728*pow(Bep,4)*Bepp*pow(Qr,2)*pow(ssr,5)*pow(P,4) + 3840*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,5)*pow(P,4) - 41856*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,5)*pow(P,4) + 72192*pow(Bep,4)*pow(ssr,7)*pow(P,4) - 724992*pow(Bep,5)*Qr*pow(ssr,7)*pow(P,4) + 35328*pow(Bep,3)*Bepp*Qr*pow(r,2)*pow(ssr,7)*pow(P,4) - 466944*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,7)*pow(P,4) - 884736*pow(Bep,4)*pow(Bepp,2)*pow(Qr,2)*pow(r,2)*pow(ssr,9)*pow(P,4) + 3200*pow(Bep,3)*pow(ssr,4)*pow(P,5) - 35200*pow(Bep,4)*Qr*pow(ssr,4)*pow(P,5) + 30720*pow(Bep,3)*Bepp*pow(ssr,6)*pow(P,5) - 337920*pow(Bep,4)*Bepp*Qr*pow(ssr,6)*pow(P,5) - 9600*pow(Bep,4)*Qr*pow(r,2)*pow(ssr,6)*pow(P,5) - 184320*pow(Bep,5)*pow(ssr,8)*pow(P,5) - 92160*pow(Bep,4)*Bepp*Qr*pow(r,2)*pow(ssr,8)*pow(P,5) - 7680*pow(Bep,4)*pow(ssr,5)*pow(P,6) - 73728*pow(Bep,4)*Bepp*pow(ssr,7)*pow(P,6)))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*pow(-1 + 8*Bep*Qr + 12*Bep*ssr*P,2)) + r_Der_P*((-12288*pow(r,5)*r_Der_Q*pow(nn,2)*(-(pow(Bep,4)*pow(ssr,6)) + 8*pow(Bep,5)*Qr*pow(ssr,6) - pow(Bep,4)*pow(r,2)*pow(ssr,8) + 8*pow(Bep,5)*pow(ssr,7)*P + 8*pow(Bep,5)*pow(r,2)*pow(ssr,9)*P))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*(-1 + 8*Bep*Qr + 12*Bep*ssr*P)) + (16*pow(r,5)*pow(nn,2)*(Bep*pow(ssr,2) - 24*pow(Bep,2)*Qr*pow(ssr,2) + 192*pow(Bep,3)*pow(Qr,2)*pow(ssr,2) - 512*pow(Bep,4)*pow(Qr,3)*pow(ssr,2) + Bep*pow(r,2)*pow(ssr,4) - 24*pow(Bep,2)*Qr*pow(r,2)*pow(ssr,4) + 224*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(ssr,4) - 768*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(ssr,4) - 96*pow(Bep,3)*pow(ssr,6) + 768*pow(Bep,4)*Qr*pow(ssr,6) + 768*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,6) - 6144*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,6) + 32*pow(Bep,3)*pow(Qr,2)*pow(r,4)*pow(ssr,6) + 128*pow(Bep,4)*pow(Qr,3)*pow(r,4)*pow(ssr,6) - 96*pow(Bep,3)*pow(r,2)*pow(ssr,8) + 768*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,4)*pow(ssr,8) - 28*pow(Bep,2)*pow(ssr,3)*P + 448*pow(Bep,3)*Qr*pow(ssr,3)*P - 1792*pow(Bep,4)*pow(Qr,2)*pow(ssr,3)*P - 28*pow(Bep,2)*pow(r,2)*pow(ssr,5)*P + 576*pow(Bep,3)*Qr*pow(r,2)*pow(ssr,5)*P - 3008*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,5)*P + 768*pow(Bep,4)*pow(ssr,7)*P + 1536*pow(Bep,3)*Bepp*Qr*pow(r,2)*pow(ssr,7)*P - 18432*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,7)*P - 192*pow(Bep,4)*pow(Qr,2)*pow(r,4)*pow(ssr,7)*P + 768*pow(Bep,4)*pow(r,2)*pow(ssr,9)*P - 6144*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,4)*pow(ssr,9)*P + 288*pow(Bep,3)*pow(ssr,4)*pow(P,2) - 2304*pow(Bep,4)*Qr*pow(ssr,4)*pow(P,2) + 288*pow(Bep,3)*pow(r,2)*pow(ssr,6)*pow(P,2) - 3456*pow(Bep,4)*Qr*pow(r,2)*pow(ssr,6)*pow(P,2) - 18432*pow(Bep,4)*Bepp*Qr*pow(r,2)*pow(ssr,8)*pow(P,2) - 960*pow(Bep,4)*pow(ssr,5)*pow(P,3) - 960*pow(Bep,4)*pow(r,2)*pow(ssr,7)*pow(P,3)))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*(-1 + 8*Bep*Qr + 12*Bep*ssr*P)));

  ans /= (3072*r_Der_Q*pow(ssr,5)*nn*(-(pow(Bep,3)*pow(r,5)) + 8*pow(Bep,4)*Qr*pow(r,5) + 8*pow(Bep,4)*pow(r,5)*ssr*P))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P) + (4*ssr*nn*(-pow(r,5) + 24*Bep*Qr*pow(r,5) - 192*pow(Bep,2)*pow(Qr,2)*pow(r,5) + 512*pow(Bep,3)*pow(Qr,3)*pow(r,5) - 32*pow(Bep,2)*pow(Qr,2)*pow(r,7)*pow(ssr,2) + 256*pow(Bep,3)*pow(Qr,3)*pow(r,7)*pow(ssr,2) + 96*pow(Bep,2)*pow(r,5)*pow(ssr,4) - 768*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4) - 768*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,7)*pow(ssr,4) + 6144*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,7)*pow(ssr,4) + 28*Bep*pow(r,5)*ssr*P - 448*pow(Bep,2)*Qr*pow(r,5)*ssr*P + 1792*pow(Bep,3)*pow(Qr,2)*pow(r,5)*ssr*P + 192*pow(Bep,3)*pow(Qr,2)*pow(r,7)*pow(ssr,3)*P - 768*pow(Bep,3)*pow(r,5)*pow(ssr,5)*P + 6144*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,7)*pow(ssr,5)*P - 288*pow(Bep,2)*pow(r,5)*pow(ssr,2)*pow(P,2) + 2304*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,2)*pow(P,2) + 960*pow(Bep,3)*pow(r,5)*pow(ssr,3)*pow(P,3)))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P);

  return ans;
}
//==============================================================================
void Evolve_scalar_field::KO_filter(Grid_data grid,
   const string type, vector<double> &rhs, const vector<double> &vec)
{ //Fourth Order Filter
   double eps= 0.3;
   double dt = grid.dt;
   int exc_i = grid.exc_i;
   int nx = grid.nx;

   for (int i=exc_i+2; i<nx-3; ++i) {
      rhs[i] -= (eps/(16*dt))*(
        vec[i+2] -4.*vec[i+1] + 6.*vec[i] -4.*vec[i-1] + vec[i-2]
      );
   }
/*---------------------------------------------------------------------------*/
   if (exc_i>0) return;
/*---------------------------------------------------------------------------*/
   if (type=="even") {
      rhs[1]-= (eps/(16*dt))*(
        vec[3] - 4.*vec[2] + 6.*vec[1] - 4.*vec[0] + vec[1]
      );
      rhs[0]-=  (eps/(16*dt))*(
        vec[2] - 4.*vec[1] + 6.*vec[0] - 4.*vec[1] + vec[2]
      );
   } else
   if (type=="odd") {
     rhs[1]-= (eps/(16*dt))*(
       vec[3] - 4.*vec[2] + 6.*vec[1] - 4.*vec[0] - vec[1]
     );
     rhs[0]-=  (eps/(16*dt))*(
       vec[2] - 4.*vec[1] + 6.*vec[0] + 4.*vec[1] - vec[2]
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
  vector<double> &dphidt)
  {

  //Field ordering p_v,q_v,phi_v
  assert(grid.exc_i ==0);
  vector<double> r = grid.r;
  int nx = grid.nx;
  vector<double> dr = grid.dr;
  double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;


  {
      int i =grid.exc_i;

      double nn = n_v.v[i];
      double r_Der_nn = 0.;
      double r_Der_nn_2 = Dx_2_ptpc_2nd(n_v.v[i+1], n_v.v[i], n_v.v[i+1], dr[i]);

      double ss = s_v.v[i];
      double r_Der_ss = Dx_ptpc_2nd(s_v.v[i+1], -s_v.v[i+1], dr[i]);
      double r_Der_ss_3 = Dx_3_ptpc_3(s_v.v[i+2], s_v.v[i+1], -s_v.v[i+1], -s_v.v[i+2], dr[i]);

      double P = p_v.v[i];
      double r_Der_P = 0.;
      double rr_Der_P =Dx_2_ptpc_2nd(p_v.v[i+1], p_v.v[i], p_v.v[i+1], dr[i]);
      double r_Der_P_2 = Dx_2_ptpc_2nd(p_v.v[i+1], p_v.v[i], p_v.v[i+1], dr[i]);

      double Q = q_v.v[i];
      double r_Der_Q = Dx_ptpc_2nd(q_v.v[i+1], -q_v.v[i+1], dr[i]);
      double r_Der_Q_3 = Dx_3_ptpc_3(q_v.v[i+2], q_v.v[i+1], -q_v.v[i+1], -q_v.v[i+2], dr[i]);


      double Bep = beta_p(ls,lexp,mu, phi_v.v[i]);
      double Bepp = beta_pp(ls,lexp,mu, phi_v.v[i]);
      double Beppp = beta_ppp(ls, lexp, mu, phi_v.v[i]);

      //Written this way because of expansion around zero and appearance of second derivatives.
      dpdt[i] = -(nn*(pow(P,3) - 16*Bep*pow(P,3)*r_Der_Q + 64*pow(Bep,2)*pow(P,3)*pow(r_Der_Q,2) - 20*Bep*pow(P,4)*r_Der_ss + 12*r_Der_Q*r_Der_ss + 160*pow(Bep,2)*pow(P,4)*r_Der_Q*r_Der_ss - 288*Bep*pow(r_Der_Q,2)*r_Der_ss + 2304*pow(Bep,2)*pow(r_Der_Q,3)*r_Der_ss - 6144*pow(Bep,3)*pow(r_Der_Q,4)*r_Der_ss + 6*P*pow(r_Der_ss,2) + 80*pow(Bep,2)*pow(P,5)*pow(r_Der_ss,2) - 480*Bep*P*r_Der_Q*pow(r_Der_ss,2) + 6528*pow(Bep,2)*P*pow(r_Der_Q,2)*pow(r_Der_ss,2) - 24576*pow(Bep,3)*P*pow(r_Der_Q,3)*pow(r_Der_ss,2) - 208*Bep*pow(P,2)*pow(r_Der_ss,3) - 128*Bep*Bepp*pow(P,4)*pow(r_Der_ss,3) + 6272*pow(Bep,2)*pow(P,2)*r_Der_Q*pow(r_Der_ss,3) + 1024*pow(Bep,2)*Bepp*pow(P,4)*r_Der_Q*pow(r_Der_ss,3) - 36864*pow(Bep,3)*pow(P,2)*pow(r_Der_Q,2)*pow(r_Der_ss,3) + 2208*pow(Bep,2)*pow(P,3)*pow(r_Der_ss,4) + 768*pow(Bep,2)*Bepp*pow(P,5)*pow(r_Der_ss,4) - 26112*pow(Bep,3)*pow(P,3)*r_Der_Q*pow(r_Der_ss,4) + 48*Bep*pow(r_Der_ss,5) + 384*Bep*Bepp*pow(P,2)*pow(r_Der_ss,5) - 7680*pow(Bep,3)*pow(P,4)*pow(r_Der_ss,5) - 768*pow(Bep,2)*r_Der_Q*pow(r_Der_ss,5) - 6144*pow(Bep,2)*Bepp*pow(P,2)*r_Der_Q*pow(r_Der_ss,5) + 3072*pow(Bep,3)*pow(r_Der_Q,2)*pow(r_Der_ss,5) + 24576*pow(Bep,3)*Bepp*pow(P,2)*pow(r_Der_Q,2)*pow(r_Der_ss,5) - 384*pow(Bep,2)*P*pow(r_Der_ss,6) - 3072*pow(Bep,2)*Bepp*pow(P,3)*pow(r_Der_ss,6) + 3072*pow(Bep,3)*P*r_Der_Q*pow(r_Der_ss,6) + 24576*pow(Bep,3)*Bepp*pow(P,3)*r_Der_Q*pow(r_Der_ss,6)));

      dpdt[i] /= 4*r_Der_ss*(-1 + 24*Bep*r_Der_Q - 192*pow(Bep,2)*pow(r_Der_Q,2) + 512*pow(Bep,3)*pow(r_Der_Q,3) + 28*Bep*P*r_Der_ss - 448*pow(Bep,2)*P*r_Der_Q*r_Der_ss + 1792*pow(Bep,3)*P*pow(r_Der_Q,2)*r_Der_ss - 288*pow(Bep,2)*pow(P,2)*pow(r_Der_ss,2) + 2304*pow(Bep,3)*pow(P,2)*r_Der_Q*pow(r_Der_ss,2) + 960*pow(Bep,3)*pow(P,3)*pow(r_Der_ss,3) + 96*pow(Bep,2)*pow(r_Der_ss,4) - 1536*pow(Bep,3)*r_Der_Q*pow(r_Der_ss,4) + 6144*pow(Bep,4)*pow(r_Der_Q,2)*pow(r_Der_ss,4) - 768*pow(Bep,3)*P*pow(r_Der_ss,5) + 6144*pow(Bep,4)*P*r_Der_Q*pow(r_Der_ss,5));

      // double p2num = 1179648*pow(Bep,5)*pow(P,2)*(-(Bep*P*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(1 + 8*Bepp*pow(P,2) - 40*pow(Bep,2)*pow(P,4) - 16*Bep*(1 + 8*Bepp*pow(P,2))*r_Der_Q + 64*pow(Bep,2)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,2))) + nn*(176*pow(Bep,3)*P*pow(r_Der_P_2,2)*pow(1 - 8*Bep*r_Der_Q,2) - Bepp*P*pow(-1 + 8*Bep*r_Der_Q,3)*(-r_Der_Q + 6*Bep*pow(r_Der_Q,2)) + Bep*r_Der_P_2*(-1 + 8*Bep*r_Der_Q)*(-3 - 16*Bepp*pow(P,2) + 200*pow(Bep,2)*pow(P,4) + 6*Bep*(11 - 16*Bepp*pow(P,2))*r_Der_Q + 32*pow(Bep,2)*(-15 + 56*Bepp*pow(P,2))*pow(r_Der_Q,2) + 1152*pow(Bep,3)*pow(r_Der_Q,3)) + 8*pow(P,3)*(-1 + 8*Bep*r_Der_Q)*((pow(Bepp,2) - Bep*Beppp)*r_Der_Q + 16*Bep*(-pow(Bepp,2) + Bep*Beppp)*pow(r_Der_Q,2) + 176*pow(Bep,2)*pow(Bepp,2)*pow(r_Der_Q,3) - 2*Bep*pow(r_Der_Q,2)*(11*pow(Bepp,2) - 32*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(1 + 32*Beppp*r_Der_Q))) - 320*pow(Bep,3)*pow(P,5)*(3*Bepp*pow(r_Der_Q,2) + Bep*r_Der_Q_3)))*pow(r_Der_ss,11) - 2359296*pow(Bep,6)*pow(P,3)*(-1 + 8*Bep*r_Der_Q)*(Bep*P*(1 + 8*Bepp*pow(P,2))*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q) - nn*(192*pow(Bep,3)*P*pow(r_Der_P_2,2) + 8*Bep*r_Der_P_2*(1 + 2*Bepp*pow(P,2) + Bep*(-15 + 32*Bepp*pow(P,2))*r_Der_Q + 56*pow(Bep,2)*pow(r_Der_Q,2)) + P*(8*Bep*Bepp*pow(r_Der_Q,2)*(1 + 24*Bepp*pow(P,2) - 8*Bep*r_Der_Q) + (Bepp + 8*pow(Bepp,2)*pow(P,2) - 8*Bep*Beppp*pow(P,2))*r_Der_Q*(-1 + 8*Bep*r_Der_Q))))*pow(r_Der_ss,12) + 37748736*pow(Bep,8)*nn*pow(P,4)*r_Der_P_2*pow(1 - 8*Bep*r_Der_Q,2)*pow(r_Der_ss,13) + (nn*pow(P,3)*pow(-1 + 8*Bep*r_Der_Q,7)*r_Der_ss_3)/(24.*pow(r_Der_ss,2)) + (P*pow(1 - 8*Bep*r_Der_Q,6)*(3*pow(P,2)*r_Der_nn_2*(1 - 8*Bep*r_Der_Q) + nn*(9*P*r_Der_P_2*(1 - 8*Bep*r_Der_Q) + 6*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + 8*pow(P,2)*(3*Bepp*pow(r_Der_Q,2) + Bep*r_Der_Q_3) + 72*Bep*pow(P,3)*r_Der_ss_3)))/(24.*r_Der_ss) + 8192*pow(Bep,4)*P*pow(r_Der_ss,10)*(-9*Bep*P*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(-3 - 24*Bepp*pow(P,2) + 464*pow(Bep,2)*pow(P,4) + 384*pow(Bep,2)*Bepp*pow(P,6) - 8*Bep*(-9 - 72*Bepp*pow(P,2) + 512*pow(Bep,2)*pow(P,4))*r_Der_Q - 576*pow(Bep,2)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,2) + 1536*pow(Bep,3)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,3)) + nn*(-27*Bepp*P*pow(1 - 8*Bep*r_Der_Q,4)*(-r_Der_Q + 4*Bep*pow(r_Der_Q,2)) + 4320*pow(Bep,3)*P*pow(r_Der_P_2,2)*(-1 + 16*pow(Bep,2)*pow(P,4) + 24*Bep*r_Der_Q - 192*pow(Bep,2)*pow(r_Der_Q,2) + 512*pow(Bep,3)*pow(r_Der_Q,3)) + 36*Bep*r_Der_P_2*(-1 - 12*Bepp*pow(P,2) + 488*pow(Bep,2)*pow(P,4) + 64*pow(Bep,2)*Bepp*pow(P,6) + Bep*(37 + 144*Bepp*pow(P,2) - 8512*pow(Bep,2)*pow(P,4) + 3328*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q + 32*pow(Bep,2)*(-17 + 36*Bepp*pow(P,2) + 1152*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,2) - 128*pow(Bep,3)*(-31 + 168*Bepp*pow(P,2))*pow(r_Der_Q,3) + 2048*pow(Bep,4)*(-7 + 36*Bepp*pow(P,2))*pow(r_Der_Q,4) + 20480*pow(Bep,5)*pow(r_Der_Q,5)) + 24*pow(P,3)*pow(1 - 8*Bep*r_Der_Q,2)*(9*(pow(Bepp,2) - Bep*Beppp)*r_Der_Q + 144*Bep*(-pow(Bepp,2) + Bep*Beppp)*pow(r_Der_Q,2) + 1440*pow(Bep,2)*pow(Bepp,2)*pow(r_Der_Q,3) - Bep*pow(r_Der_Q,2)*(180*pow(Bepp,2) - 576*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(71 + 576*Beppp*r_Der_Q))) - 1152*pow(Bep,2)*pow(P,7)*(-28*Bep*pow(Bepp,2)*pow(r_Der_Q,2) + r_Der_Q*(pow(Bepp,2)*(3 - 8*Bep*r_Der_Q) + 3*Bep*Beppp*(-1 + 8*Bep*r_Der_Q)) + 16*pow(Bep,2)*Bepp*r_Der_Q_3) - 16*pow(Bep,2)*pow(P,5)*(3*Bepp*(16*Bep*pow(r_Der_Q,2)*(-67 + 728*Bep*r_Der_Q) + 3*r_Der_Q*(9 - 264*Bep*r_Der_Q + 2048*pow(Bep,2)*pow(r_Der_Q,2))) + 128*pow(Bep,2)*(-11 + 124*Bep*r_Der_Q)*r_Der_Q_3) - 1920*pow(Bep,4)*pow(P,6)*(-1 + 8*Bep*r_Der_Q)*r_Der_ss_3)) + (8*Bep*pow(1 - 8*Bep*r_Der_Q,2)*pow(r_Der_ss,3)*(3*Bep*pow(P,3)*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(-251 - 100*Bepp*pow(P,2) - 3828*pow(Bep,2)*pow(P,4) + 320*Bep*(43 + 5*Bepp*pow(P,2))*r_Der_Q - 64*pow(Bep,2)*(2687 + 100*Bepp*pow(P,2))*pow(r_Der_Q,2) + 623616*pow(Bep,3)*pow(r_Der_Q,3)) + nn*(P*(-26604*pow(Bep,3)*pow(P,5)*r_Der_P_2*(-1 + 8*Bep*r_Der_Q) - 1068*Bep*Bepp*pow(P,3)*r_Der_P_2*pow(-1 + 8*Bep*r_Der_Q,3) - 3*Bep*pow(r_Der_Q,2)*pow(-1 + 8*Bep*r_Der_Q,3)*(-51 + 1208*Bep*r_Der_Q) + 3*Bep*P*r_Der_P_2*pow(-1 + 8*Bep*r_Der_Q,3)*(-563 + 5632*Bep*r_Der_Q) + 24*pow(Bep,2)*pow(P,6)*(3*Bepp*r_Der_Q*(-109 + 1276*Bep*r_Der_Q) + 404*pow(Bep,2)*r_Der_Q_3) + 4*pow(P,2)*pow(1 - 8*Bep*r_Der_Q,2)*(3*Bepp*(146*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + r_Der_Q*(7 + 94*Bep*r_Der_Q + 48*pow(Bep,2)*pow(r_Der_Q,2))) + 4*pow(Bep,2)*(-737 + 6208*Bep*r_Der_Q)*r_Der_Q_3) - 12*pow(P,4)*(-1 + 8*Bep*r_Der_Q)*((14*pow(Bepp,2) + 25*Bep*Beppp)*r_Der_Q + 2*pow(Bep,2)*pow(r_Der_Q,2)*(-800*pow(Bepp,2)*r_Der_Q + Bep*(-1613 + 800*Beppp*r_Der_Q)) + 104*pow(Bep,2)*Bepp*r_Der_Q_3 - 8*r_Der_Q*(Bep*(-11*pow(Bepp,2) + 50*Bep*Beppp)*r_Der_Q + 104*pow(Bep,3)*Bepp*r_Der_Q_3))) - 3*(-1 - 8*Bepp*pow(P,2) + 1880*pow(Bep,2)*pow(P,4) + 800*pow(Bep,2)*Bepp*pow(P,6) - 12640*pow(Bep,4)*pow(P,8) - 16*Bep*(-1 - 20*Bepp*pow(P,2) + 2620*pow(Bep,2)*pow(P,4) + 800*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q + 320*pow(Bep,2)*(1 - 16*Bepp*pow(P,2) + 968*pow(Bep,2)*pow(P,4) + 160*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) - 10240*pow(Bep,3)*(1 - 4*Bepp*pow(P,2) + 74*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) - 20480*pow(Bep,4)*(-5 + 8*Bepp*pow(P,2))*pow(r_Der_Q,4) + 65536*pow(Bep,5)*(-7 + 4*Bepp*pow(P,2))*pow(r_Der_Q,5) + 786432*pow(Bep,6)*pow(r_Der_Q,6))*r_Der_ss_3)))/3. + pow(r_Der_ss,2)*(Bep*pow(P,2)*r_Der_nn_2*pow(1 - 8*Bep*r_Der_Q,4)*(-59 - 16*Bepp*pow(P,2) - 3360*pow(Bep,2)*pow(P,4) + 32*Bep*(151 + 8*Bepp*pow(P,2))*r_Der_Q - 64*pow(Bep,2)*(1031 + 16*Bepp*pow(P,2))*pow(r_Der_Q,2) + 248832*pow(Bep,3)*pow(r_Der_Q,3)) + (nn*pow(-1 + 8*Bep*r_Der_Q,3)*(-6*Bep*pow(r_Der_Q,2)*pow(-1 + 8*Bep*r_Der_Q,3)*(-3 + 88*Bep*r_Der_Q) + 24*Bep*P*r_Der_P_2*pow(-1 + 8*Bep*r_Der_Q,3)*(-17 + 168*Bep*r_Der_Q) + 32*pow(Bep,2)*pow(P,6)*(3*Bepp*r_Der_Q*(-53 + 840*Bep*r_Der_Q) + 416*pow(Bep,2)*r_Der_Q_3) + pow(P,2)*pow(1 - 8*Bep*r_Der_Q,2)*(3*Bepp*(80*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + r_Der_Q*(5 + 88*Bep*r_Der_Q)) + 16*pow(Bep,2)*(-337 + 2760*Bep*r_Der_Q)*r_Der_Q_3) - 16*pow(P,4)*(-1 + 8*Bep*r_Der_Q)*(3*(pow(Bepp,2) + Bep*Beppp)*r_Der_Q + 3*pow(Bep,2)*pow(r_Der_Q,2)*(-64*pow(Bepp,2)*r_Der_Q + Bep*(-617 + 64*Beppp*r_Der_Q)) + 16*pow(Bep,2)*Bepp*r_Der_Q_3 - 16*pow(Bep,2)*r_Der_Q*(3*Beppp*r_Der_Q + 8*Bep*Bepp*r_Der_Q_3)) + 51552*pow(Bep,4)*pow(P,7)*r_Der_ss_3 - 96*pow(Bep,2)*pow(P,5)*(-1 + 8*Bep*r_Der_Q)*(262*Bep*r_Der_P_2 + 13*Bepp*(-1 + 8*Bep*r_Der_Q)*r_Der_ss_3) + 24*Bep*pow(P,3)*pow(1 - 8*Bep*r_Der_Q,2)*(8*Bepp*r_Der_P_2*(1 - 8*Bep*r_Der_Q) + Bep*(-95 + 552*Bep*r_Der_Q)*r_Der_ss_3)))/3.) - 64*pow(Bep,2)*pow(r_Der_ss,8)*(3*Bep*r_Der_nn_2*(3 + 24*Bepp*pow(P,2) - 9376*pow(Bep,2)*pow(P,4) - 25856*pow(Bep,2)*Bepp*pow(P,6) + 254720*pow(Bep,4)*pow(P,8) + 30720*pow(Bep,4)*Bepp*pow(P,10) - 16*Bep*(9 + 72*Bepp*pow(P,2) - 20976*pow(Bep,2)*pow(P,4) - 38784*pow(Bep,2)*Bepp*pow(P,6) + 148480*pow(Bep,4)*pow(P,8))*r_Der_Q - 192*pow(Bep,2)*(-15 - 120*Bepp*pow(P,2) + 23200*pow(Bep,2)*pow(P,4) + 25856*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 2048*pow(Bep,3)*(-15 - 120*Bepp*pow(P,2) + 12712*pow(Bep,2)*pow(P,4) + 6464*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,3) - 36864*pow(Bep,4)*(-5 - 40*Bepp*pow(P,2) + 1536*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,4) - 589824*pow(Bep,5)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,5) + 786432*pow(Bep,6)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,6)) + nn*(92160*pow(Bep,4)*(-pow(Bepp,2) + Bep*Beppp)*pow(P,10)*r_Der_Q - 144*Bep*Bepp*P*r_Der_P_2*pow(-1 + 8*Bep*r_Der_Q,5)*(1 + 8*Bep*r_Der_Q) - 192*pow(Bep,3)*pow(P,3)*r_Der_P_2*pow(-1 + 8*Bep*r_Der_Q,3)*(-369 + 5152*Bep*r_Der_Q) - 9*pow(-1 + 8*Bep*r_Der_Q,5)*(128*pow(Bep,3)*pow(r_Der_P_2,2) + Bepp*r_Der_Q*(-1 + 8*Bep*r_Der_Q)) - 24*pow(P,2)*pow(1 - 8*Bep*r_Der_Q,4)*(3*(pow(Bepp,2) - Bep*Beppp)*r_Der_Q + 48*Bep*(-pow(Bepp,2) + Bep*Beppp)*pow(r_Der_Q,2) + 384*pow(Bep,2)*pow(Bepp,2)*pow(r_Der_Q,3) - 2*Bep*pow(r_Der_Q,2)*(24*pow(Bepp,2) - 96*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(179 + 96*Beppp*r_Der_Q))) - 256*pow(Bep,4)*pow(P,8)*(3*Bepp*(1368*Bep*pow(r_Der_Q,2) + 5*r_Der_Q*(-5 + 32*Bep*r_Der_Q)) + 1760*pow(Bep,2)*r_Der_Q_3) + 96*pow(Bep,2)*pow(P,4)*pow(1 - 8*Bep*r_Der_Q,2)*(-15264*pow(Bep,3)*pow(r_Der_P_2,2) + Bepp*(16*Bep*pow(r_Der_Q,2)*(-93 + 3304*Bep*r_Der_Q) + r_Der_Q*(209 - 4104*Bep*r_Der_Q + 27648*pow(Bep,2)*pow(r_Der_Q,2))) + 32*pow(Bep,2)*(9 + 440*Bep*r_Der_Q)*r_Der_Q_3) + 768*pow(Bep,2)*pow(P,6)*((53*pow(Bepp,2) - 101*Bep*Beppp)*r_Der_Q - 40192*pow(Bep,3)*pow(Bepp,2)*pow(r_Der_Q,4) + 8*pow(Bep,2)*pow(r_Der_Q,3)*(1256*pow(Bepp,2) - 1344*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(797 + 6464*Beppp*r_Der_Q)) + 512*pow(Bep,2)*Bepp*r_Der_Q_3 - 8*r_Der_Q*(Bep*(127*pow(Bepp,2) - 303*Bep*Beppp)*r_Der_Q + 1024*pow(Bep,3)*Bepp*r_Der_Q_3) + Bep*pow(r_Der_Q,2)*(-628*pow(Bepp,2) + 6080*Bep*pow(Bepp,2)*r_Der_Q - pow(Bep,2)*(701 + 19392*Beppp*r_Der_Q) + 32768*pow(Bep,3)*Bepp*r_Der_Q_3)) - 61440*pow(Bep,5)*pow(P,9)*(-3*Bepp*r_Der_P_2 + 5*Bep*r_Der_ss_3) - 6144*pow(Bep,4)*pow(P,7)*(5*Bep*r_Der_P_2*(-50 + 451*Bep*r_Der_Q) - 26*Bepp*pow(1 - 8*Bep*r_Der_Q,2)*r_Der_ss_3) + 1536*pow(Bep,3)*pow(P,5)*pow(1 - 8*Bep*r_Der_Q,2)*(-9*Bepp*r_Der_P_2*(5 + 172*Bep*r_Der_Q) + Bep*(-29 + 648*Bep*r_Der_Q)*r_Der_ss_3))) - 256*pow(Bep,3)*P*pow(r_Der_ss,7)*(-6*Bep*pow(P,2)*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(119 + 460*Bepp*pow(P,2) - 11512*pow(Bep,2)*pow(P,4) - 2752*pow(Bep,2)*Bepp*pow(P,6) + 4*Bep*(-1137 - 2760*Bepp*pow(P,2) + 32512*pow(Bep,2)*pow(P,4))*r_Der_Q + 96*pow(Bep,2)*(661 + 920*Bepp*pow(P,2))*pow(r_Der_Q,2) - 256*pow(Bep,3)*(1507 + 920*Bepp*pow(P,2))*pow(r_Der_Q,3) + 866304*pow(Bep,4)*pow(r_Der_Q,4)) + nn*(144*Bep*pow(r_Der_Q,2)*pow(-1 + 8*Bep*r_Der_Q,5) - 9*Bep*P*r_Der_P_2*pow(1 - 8*Bep*r_Der_Q,4)*(-135 + 2408*Bep*r_Der_Q) - 128*pow(Bep,2)*pow(P,8)*(3*r_Der_Q*(43*Bep*Beppp*(1 - 8*Bep*r_Der_Q) + pow(Bepp,2)*(-35 + 344*Bep*r_Der_Q)) + 64*pow(Bep,2)*Bepp*r_Der_Q_3) - 16*pow(Bep,2)*pow(P,6)*(3*Bepp*(3484*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + r_Der_Q*(65 - 664*Bep*r_Der_Q + 3200*pow(Bep,2)*pow(r_Der_Q,2))) + 128*pow(Bep,2)*(-51 + 424*Bep*r_Der_Q)*r_Der_Q_3) + pow(P,2)*pow(-1 + 8*Bep*r_Der_Q,3)*(-50880*pow(Bep,3)*pow(r_Der_P_2,2) + 3*Bepp*(32*Bep*pow(r_Der_Q,2)*(-9 + 1672*Bep*r_Der_Q) + r_Der_Q*(167 - 4464*Bep*r_Der_Q + 27072*pow(Bep,2)*pow(r_Der_Q,2))) + 8*pow(Bep,2)*(427 + 3240*Bep*r_Der_Q)*r_Der_Q_3) - 8*pow(P,4)*(-1 + 8*Bep*r_Der_Q)*(-69*(pow(Bepp,2) - 5*Bep*Beppp)*r_Der_Q + 99840*pow(Bep,3)*pow(Bepp,2)*pow(r_Der_Q,4) - 24*pow(Bep,2)*pow(r_Der_Q,3)*(1040*pow(Bepp,2) - 960*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(2327 + 7360*Beppp*r_Der_Q)) - 1664*pow(Bep,2)*Bepp*r_Der_Q_3 + Bep*pow(r_Der_Q,2)*(1560*pow(Bepp,2) - 10176*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(5253 + 66240*Beppp*r_Der_Q) - 106496*pow(Bep,3)*Bepp*r_Der_Q_3) + 8*r_Der_Q*(3*Bep*(61*pow(Bepp,2) - 345*Bep*Beppp)*r_Der_Q + 3328*pow(Bep,3)*Bepp*r_Der_Q_3)) - 256*pow(Bep,3)*pow(P,7)*(-1 + 8*Bep*r_Der_Q)*(-141*Bepp*r_Der_P_2 + 265*Bep*r_Der_ss_3) + 16*pow(Bep,2)*pow(P,5)*(-1 + 8*Bep*r_Der_Q)*(3*Bep*r_Der_P_2*(2911 - 28580*Bep*r_Der_Q) + 832*Bepp*pow(1 - 8*Bep*r_Der_Q,2)*r_Der_ss_3) + 16*Bep*pow(P,3)*pow(-1 + 8*Bep*r_Der_Q,3)*(-3*Bepp*r_Der_P_2*(61 + 1632*Bep*r_Der_Q) + 2*Bep*(-35 + 1944*Bep*r_Der_Q)*r_Der_ss_3))) + (64*pow(Bep,2)*pow(r_Der_ss,6)*(-18*Bep*pow(P,2)*r_Der_nn_2*(-29 - 148*Bepp*pow(P,2) + 9502*pow(Bep,2)*pow(P,4) + 3248*pow(Bep,2)*Bepp*pow(P,6) + 1600*pow(Bep,4)*pow(P,8) - 32*Bep*(-52 - 185*Bepp*pow(P,2) + 8981*pow(Bep,2)*pow(P,4) + 1624*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q + 64*pow(Bep,2)*(-605 - 1480*Bepp*pow(P,2) + 43342*pow(Bep,2)*pow(P,4) + 3248*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) - 20480*pow(Bep,3)*(-23 - 37*Bepp*pow(P,2) + 423*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) - 20480*pow(Bep,4)*(155 + 148*Bepp*pow(P,2))*pow(r_Der_Q,4) + 131072*pow(Bep,5)*(86 + 37*Bepp*pow(P,2))*pow(r_Der_Q,5) - 16515072*pow(Bep,6)*pow(r_Der_Q,6)) + nn*(28800*pow(Bep,4)*Bepp*pow(P,10)*r_Der_Q - 45*Bep*pow(r_Der_Q,2)*pow(1 - 8*Bep*r_Der_Q,6) + 504*Bep*P*r_Der_P_2*pow(-1 + 8*Bep*r_Der_Q,5)*(-1 + 24*Bep*r_Der_Q) - 6*pow(P,2)*pow(1 - 8*Bep*r_Der_Q,4)*(-7968*pow(Bep,3)*pow(r_Der_P_2,2) + 3*Bepp*(32*Bep*pow(r_Der_Q,2)*(3 + 256*Bep*r_Der_Q) + r_Der_Q*(13 - 752*Bep*r_Der_Q + 4032*pow(Bep,2)*pow(r_Der_Q,2))) + 8*pow(Bep,2)*(97 + 200*Bep*r_Der_Q)*r_Der_Q_3) + 4*pow(Bep,2)*pow(P,6)*(-1 + 8*Bep*r_Der_Q)*(3*Bepp*(43520*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + r_Der_Q*(927 - 3672*Bep*r_Der_Q + 37632*pow(Bep,2)*pow(r_Der_Q,2))) + 64*pow(Bep,2)*(-1957 + 16712*Bep*r_Der_Q)*r_Der_Q_3) - 96*pow(Bep,2)*pow(P,8)*(3*(-115*pow(Bepp,2) + 203*Bep*Beppp)*r_Der_Q + 12*pow(Bep,2)*pow(r_Der_Q,2)*(-3248*pow(Bepp,2)*r_Der_Q + Bep*(-125 + 3248*Beppp*r_Der_Q)) + 704*pow(Bep,2)*Bepp*r_Der_Q_3 - 16*r_Der_Q*(3*Bep*(-159*pow(Bepp,2) + 203*Bep*Beppp)*r_Der_Q + 352*pow(Bep,3)*Bepp*r_Der_Q_3)) + 24*pow(P,4)*pow(1 - 8*Bep*r_Der_Q,2)*(3*(7*pow(Bepp,2) + 37*Bep*Beppp)*r_Der_Q + 19968*pow(Bep,3)*pow(Bepp,2)*pow(r_Der_Q,4) - 8*pow(Bep,2)*pow(r_Der_Q,3)*(624*pow(Bepp,2) - 384*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(5779 + 7104*Beppp*r_Der_Q)) - 488*pow(Bep,2)*Bepp*r_Der_Q_3 + Bep*pow(r_Der_Q,2)*(312*pow(Bepp,2) + 576*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(3635 + 21312*Beppp*r_Der_Q) - 31232*pow(Bep,3)*Bepp*r_Der_Q_3) + 8*r_Der_Q*(-9*Bep*(4*pow(Bepp,2) + 37*Bep*Beppp)*r_Der_Q + 976*pow(Bep,3)*Bepp*r_Der_Q_3)) - 48*pow(Bep,2)*pow(P,5)*pow(1 - 8*Bep*r_Der_Q,2)*(3*Bep*r_Der_P_2*(2435 - 25088*Bep*r_Der_Q) + 584*Bepp*pow(1 - 8*Bep*r_Der_Q,2)*r_Der_ss_3) - 768*pow(Bep,4)*pow(P,9)*(75*Bep*r_Der_P_2 + 16*Bepp*(-1 + 8*Bep*r_Der_Q)*r_Der_ss_3) + 192*pow(Bep,3)*pow(P,7)*(-1 + 8*Bep*r_Der_Q)*(741*Bepp*r_Der_P_2*(1 - 8*Bep*r_Der_Q) + 4*Bep*(-401 + 3144*Bep*r_Der_Q)*r_Der_ss_3) - 12*Bep*pow(P,3)*pow(1 - 8*Bep*r_Der_Q,4)*(-96*Bepp*r_Der_P_2*(3 + 59*Bep*r_Der_Q) + Bep*(25 + 9144*Bep*r_Der_Q)*r_Der_ss_3))))/3. - (2*(-1 + 8*Bep*r_Der_Q)*pow(r_Der_ss,4)*(-3*Bep*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(3 + 24*Bepp*pow(P,2) - 19448*pow(Bep,2)*pow(P,4) - 8512*pow(Bep,2)*Bepp*pow(P,6) - 82880*pow(Bep,4)*pow(P,8) + 64*Bep*(-3 - 15*Bepp*pow(P,2) + 12767*pow(Bep,2)*pow(P,4) + 2128*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q - 64*pow(Bep,2)*(-75 - 240*Bepp*pow(P,2) + 145928*pow(Bep,2)*pow(P,4) + 8512*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 61440*pow(Bep,3)*(-1 - 2*Bepp*pow(P,2) + 527*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) + 61440*pow(Bep,4)*(7 + 8*Bepp*pow(P,2))*pow(r_Der_Q,4) - 786432*pow(Bep,5)*(2 + Bepp*pow(P,2))*pow(r_Der_Q,5) + 2359296*pow(Bep,6)*pow(r_Der_Q,6)) + nn*(-960*pow(Bep,4)*Bepp*pow(P,8)*r_Der_Q*(-215 + 2072*Bep*r_Der_Q) - 24*pow(Bep,2)*Bepp*pow(P,4)*pow(1 - 8*Bep*r_Der_Q,2)*(7216*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + r_Der_Q*(257 + 2216*Bep*r_Der_Q + 4224*pow(Bep,2)*pow(r_Der_Q,2))) + 72*pow(P,2)*pow(-1 + 8*Bep*r_Der_Q,3)*(-((pow(Bepp,2) + Bep*Beppp)*r_Der_Q) + 8*Bep*(2*pow(Bepp,2) + 3*Bep*Beppp)*pow(r_Der_Q,2) + 128*pow(Bep,4)*pow(r_Der_Q,3)*(37 + 4*Beppp*r_Der_Q) - 8*pow(Bep,2)*pow(r_Der_Q,2)*(8*pow(Bepp,2)*r_Der_Q + Bep*(31 + 24*Beppp*r_Der_Q))) + 192*pow(Bep,2)*pow(P,6)*(-1 + 8*Bep*r_Der_Q)*((17*pow(Bepp,2) + 133*Bep*Beppp)*r_Der_Q + 16*Bep*(58*pow(Bepp,2) - 133*Bep*Beppp)*pow(r_Der_Q,2) + 14*pow(Bep,2)*pow(r_Der_Q,2)*(-608*pow(Bepp,2)*r_Der_Q + Bep*(-357 + 608*Beppp*r_Der_Q))) + 9*pow(-1 + 8*Bep*r_Der_Q,5)*(-128*pow(Bep,3)*pow(r_Der_P_2,2) + Bepp*(16*Bep*pow(r_Der_Q,2)*(1 + 24*Bep*r_Der_Q) + r_Der_Q*(-1 - 48*Bep*r_Der_Q + 192*pow(Bep,2)*pow(r_Der_Q,2)))) - 16*pow(Bep,2)*(9 + 12*Bepp*pow(P,2) - 32324*pow(Bep,2)*pow(P,4) + 4800*pow(Bep,2)*Bepp*pow(P,6) + 7040*pow(Bep,4)*pow(P,8) - 96*Bep*(4 + 5*Bepp*pow(P,2) - 8281*pow(Bep,2)*pow(P,4) + 800*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q + 192*pow(Bep,2)*(35 + 40*Bepp*pow(P,2) - 33924*pow(Bep,2)*pow(P,4) + 1600*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 2048*pow(Bep,3)*(-30 - 30*Bepp*pow(P,2) + 8681*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) + 61440*pow(Bep,4)*(5 + 4*Bepp*pow(P,2))*pow(r_Der_Q,4) - 393216*pow(Bep,5)*(2 + Bepp*pow(P,2))*pow(r_Der_Q,5) + 786432*pow(Bep,6)*pow(r_Der_Q,6))*r_Der_Q_3 - 476160*pow(Bep,6)*pow(P,9)*r_Der_ss_3 + 96*pow(Bep,2)*pow(P,3)*pow(-1 + 8*Bep*r_Der_Q,3)*(-3*Bep*r_Der_P_2*(-439 + 4512*Bep*r_Der_Q) + 47*Bepp*pow(1 - 8*Bep*r_Der_Q,2)*r_Der_ss_3) + 128*pow(Bep,4)*pow(P,7)*(-1 + 8*Bep*r_Der_Q)*(4215*Bep*r_Der_P_2 + 568*Bepp*(-1 + 8*Bep*r_Der_Q)*r_Der_ss_3) + 24*Bep*P*pow(-1 + 8*Bep*r_Der_Q,5)*(-6*Bepp*r_Der_P_2*(1 + 8*Bep*r_Der_Q) + Bep*(19 + 600*Bep*r_Der_Q)*r_Der_ss_3) - 32*pow(Bep,3)*pow(P,5)*pow(1 - 8*Bep*r_Der_Q,2)*(-2496*Bepp*r_Der_P_2*(-1 + 8*Bep*r_Der_Q) + Bep*(-8555 + 59352*Bep*r_Der_Q)*r_Der_ss_3))))/3. + 1024*pow(Bep,3)*pow(r_Der_ss,9)*(6*Bep*P*r_Der_nn_2*(3 + 24*Bepp*pow(P,2) - 1720*pow(Bep,2)*pow(P,4) - 3008*pow(Bep,2)*Bepp*pow(P,6) + 9600*pow(Bep,4)*pow(P,8) + 8*Bep*(-15 - 120*Bepp*pow(P,2) + 5600*pow(Bep,2)*pow(P,4) + 6016*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q - 128*pow(Bep,2)*(-15 - 120*Bepp*pow(P,2) + 3020*pow(Bep,2)*pow(P,4) + 1504*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 15360*pow(Bep,3)*(-1 - 8*Bepp*pow(P,2) + 72*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) + 61440*pow(Bep,4)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,4) - 98304*pow(Bep,5)*(1 + 8*Bepp*pow(P,2))*pow(r_Der_Q,5)) + nn*(288*pow(Bep,3)*P*pow(r_Der_P_2,2)*(-1 + 8*Bep*r_Der_Q)*(-9 + 1216*pow(Bep,2)*pow(P,4) + 216*Bep*r_Der_Q - 1728*pow(Bep,2)*pow(r_Der_Q,2) + 4608*pow(Bep,3)*pow(r_Der_Q,3)) + 3*Bep*r_Der_P_2*(3 + 96*Bepp*pow(P,2) - 11584*pow(Bep,2)*pow(P,4) - 4608*pow(Bep,2)*Bepp*pow(P,6) + 38400*pow(Bep,4)*pow(P,8) - 4*Bep*(33 + 528*Bepp*pow(P,2) - 79168*pow(Bep,2)*pow(P,4) + 39936*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q + 32*pow(Bep,2)*(75 + 192*Bepp*pow(P,2) - 88832*pow(Bep,2)*pow(P,4) + 49152*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 1536*pow(Bep,3)*(-15 + 112*Bepp*pow(P,2) + 5472*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) - 24576*pow(Bep,4)*(-5 + 64*Bepp*pow(P,2))*pow(r_Der_Q,4) + 49152*pow(Bep,5)*(-7 + 80*Bepp*pow(P,2))*pow(r_Der_Q,5) + 393216*pow(Bep,6)*pow(r_Der_Q,6)) - 2*P*(9*Bepp*pow(-1 + 8*Bep*r_Der_Q,5)*(-r_Der_Q + 2*Bep*pow(r_Der_Q,2)) - 72*pow(P,2)*pow(-1 + 8*Bep*r_Der_Q,3)*((pow(Bepp,2) - Bep*Beppp)*r_Der_Q + 16*Bep*(-pow(Bepp,2) + Bep*Beppp)*pow(r_Der_Q,2) + 144*pow(Bep,2)*pow(Bepp,2)*pow(r_Der_Q,3) - Bep*pow(r_Der_Q,2)*(18*pow(Bepp,2) - 64*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(27 + 64*Beppp*r_Der_Q))) + 8*pow(Bep,2)*pow(P,4)*(-1 + 8*Bep*r_Der_Q)*(3*Bepp*(16*Bep*pow(r_Der_Q,2)*(-117 + 2024*Bep*r_Der_Q) + r_Der_Q*(121 - 2360*Bep*r_Der_Q + 17280*pow(Bep,2)*pow(r_Der_Q,2))) + 16*pow(Bep,2)*(-97 + 2248*Bep*r_Der_Q)*r_Der_Q_3) + 64*pow(Bep,2)*pow(P,6)*(3*(-37*pow(Bepp,2) + 47*Bep*Beppp)*r_Der_Q - 8832*pow(Bep,2)*pow(Bepp,2)*pow(r_Der_Q,3) + 24*Bep*pow(r_Der_Q,2)*(46*pow(Bepp,2) - 104*Bep*pow(Bepp,2)*r_Der_Q + pow(Bep,2)*(15 + 376*Beppp*r_Der_Q)) - 736*pow(Bep,2)*Bepp*r_Der_Q_3 + 16*r_Der_Q*(3*Bep*(25*pow(Bepp,2) - 47*Bep*Beppp)*r_Der_Q + 368*pow(Bep,3)*Bepp*r_Der_Q_3)) + 6144*pow(Bep,4)*Bepp*pow(P,7)*(-1 + 8*Bep*r_Der_Q)*r_Der_ss_3 + 768*pow(Bep,4)*pow(P,5)*(7 - 144*Bep*r_Der_Q + 704*pow(Bep,2)*pow(r_Der_Q,2))*r_Der_ss_3))) + (pow(-1 + 8*Bep*r_Der_Q,5)*(6*r_Der_nn_2*(16*Bep*pow(P,4) - (3 + 128*pow(Bep,2)*pow(P,4))*r_Der_Q + 72*Bep*pow(r_Der_Q,2) - 576*pow(Bep,2)*pow(r_Der_Q,3) + 1536*pow(Bep,3)*pow(r_Der_Q,4)) + nn*(2*(-5 + 112*pow(Bep,2)*pow(P,4) + 120*Bep*r_Der_Q - 960*pow(Bep,2)*pow(r_Der_Q,2) + 2560*pow(Bep,3)*pow(r_Der_Q,3))*r_Der_Q_3 + 3*P*(4*P*(-23*Bep*P*r_Der_P_2*(-1 + 8*Bep*r_Der_Q) + 18*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + Bepp*pow(P,2)*r_Der_Q*(-1 + 64*Bep*r_Der_Q)) + (-1 + 384*pow(Bep,2)*pow(P,4) + 24*Bep*r_Der_Q - 192*pow(Bep,2)*pow(r_Der_Q,2) + 512*pow(Bep,3)*pow(r_Der_Q,3))*r_Der_ss_3))))/12. + (16*Bep*pow(r_Der_ss,5)*(12*Bep*P*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(6 + 39*Bepp*pow(P,2) - 7355*pow(Bep,2)*pow(P,4) - 3040*pow(Bep,2)*Bepp*pow(P,6) - 7680*pow(Bep,4)*pow(P,8) + 2*Bep*(-183 - 780*Bepp*pow(P,2) + 127828*pow(Bep,2)*pow(P,4) + 24320*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q - 16*pow(Bep,2)*(-555 - 1560*Bepp*pow(P,2) + 167396*pow(Bep,2)*pow(P,4) + 12160*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 768*pow(Bep,3)*(-145 - 260*Bepp*pow(P,2) + 11498*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) + 30720*pow(Bep,4)*(25 + 26*Bepp*pow(P,2))*pow(r_Der_Q,4) - 24576*pow(Bep,5)*(113 + 52*Bepp*pow(P,2))*pow(r_Der_Q,5) + 4128768*pow(Bep,6)*pow(r_Der_Q,6)) + nn*(8064*pow(Bep,3)*P*pow(r_Der_P_2,2)*pow(-1 + 8*Bep*r_Der_Q,5) + 3*Bep*r_Der_P_2*(-1 + 8*Bep*r_Der_Q)*(-9*(-1 - 28*Bepp*pow(P,2) + 6860*pow(Bep,2)*pow(P,4) + 3712*pow(Bep,2)*Bepp*pow(P,6) + 7040*pow(Bep,4)*pow(P,8)) + 32*Bep*(-21 - 147*Bepp*pow(P,2) + 51005*pow(Bep,2)*pow(P,4) + 16704*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q - 64*pow(Bep,2)*(-285 + 168*Bepp*pow(P,2) + 222820*pow(Bep,2)*pow(P,4) + 33408*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 2048*pow(Bep,3)*(-120 + 378*Bepp*pow(P,2) + 20135*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) - 12288*pow(Bep,4)*(-145 + 476*Bepp*pow(P,2))*pow(r_Der_Q,4) + 393216*pow(Bep,5)*(-17 + 35*Bepp*pow(P,2))*pow(r_Der_Q,5) + 10223616*pow(Bep,6)*pow(r_Der_Q,6)) + P*(8*pow(Bep,2)*(117 + 216*Bepp*pow(P,2) - 56896*pow(Bep,2)*pow(P,4) + 9088*pow(Bep,2)*Bepp*pow(P,6) + 1920*pow(Bep,4)*pow(P,8) - 16*Bep*(297 + 540*Bepp*pow(P,2) - 87616*pow(Bep,2)*pow(P,4) + 9088*pow(Bep,2)*Bepp*pow(P,6))*r_Der_Q + 64*pow(Bep,2)*(1215 + 2160*Bepp*pow(P,2) - 179776*pow(Bep,2)*pow(P,4) + 9088*pow(Bep,2)*Bepp*pow(P,6))*pow(r_Der_Q,2) + 30720*pow(Bep,3)*(-21 - 36*Bepp*pow(P,2) + 1024*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) + 552960*pow(Bep,4)*(5 + 8*Bepp*pow(P,2))*pow(r_Der_Q,4) - 1769472*pow(Bep,5)*(3 + 4*Bepp*pow(P,2))*pow(r_Der_Q,5) + 2359296*pow(Bep,6)*pow(r_Der_Q,6))*r_Der_Q_3 + 3*(1920*pow(Bep,4)*Bepp*pow(P,8)*r_Der_Q*(-15 + 128*Bep*r_Der_Q) - 3*Bepp*pow(-1 + 8*Bep*r_Der_Q,5)*(8*Bep*pow(r_Der_Q,2)*(9 + 344*Bep*r_Der_Q) + r_Der_Q*(-1 - 288*Bep*r_Der_Q + 1344*pow(Bep,2)*pow(r_Der_Q,2))) + 16*pow(Bep,2)*Bepp*pow(P,4)*pow(1 - 8*Bep*r_Der_Q,2)*(5940*Bep*pow(r_Der_Q,2)*(-1 + 8*Bep*r_Der_Q) + r_Der_Q*(159 + 436*Bep*r_Der_Q + 4512*pow(Bep,2)*pow(r_Der_Q,2))) + 4*pow(P,2)*pow(-1 + 8*Bep*r_Der_Q,3)*(3*(8*pow(Bepp,2) + 13*Bep*Beppp)*r_Der_Q - 24*Bep*(16*pow(Bepp,2) + 39*Bep*Beppp)*pow(r_Der_Q,2) + 3072*pow(Bep,3)*pow(Bepp,2)*pow(r_Der_Q,4) - 256*pow(r_Der_Q,3)*(3*pow(Bep,2)*pow(Bepp,2) + 26*pow(Bep,4)*(7 + 3*Beppp*r_Der_Q)) + 16*Bep*pow(r_Der_Q,2)*(3*pow(Bepp,2) + 96*Bep*pow(Bepp,2)*r_Der_Q + 4*pow(Bep,2)*(47 + 117*Beppp*r_Der_Q))) - 64*pow(Bep,2)*pow(P,6)*(-1 + 8*Bep*r_Der_Q)*(2*(-24*pow(Bepp,2) + 95*Bep*Beppp)*r_Der_Q + 16*Bep*(119*pow(Bepp,2) - 190*Bep*Beppp)*pow(r_Der_Q,2) + pow(Bep,2)*pow(r_Der_Q,2)*(-12160*pow(Bepp,2)*r_Der_Q + Bep*(-2113 + 12160*Beppp*r_Der_Q))) + 25600*pow(Bep,6)*pow(P,9)*r_Der_ss_3 - 11264*pow(Bep,4)*Bepp*pow(P,7)*pow(1 - 8*Bep*r_Der_Q,2)*r_Der_ss_3 - 3648*pow(Bep,2)*Bepp*pow(P,3)*pow(-1 + 8*Bep*r_Der_Q,5)*r_Der_ss_3 - 48*pow(Bep,2)*P*pow(-1 + 8*Bep*r_Der_Q,5)*(5 + 264*Bep*r_Der_Q)*r_Der_ss_3 + 64*pow(Bep,4)*pow(P,5)*pow(1 - 8*Bep*r_Der_Q,2)*(-1363 + 10200*Bep*r_Der_Q)*r_Der_ss_3)))))/3. + (pow(1 - 8*Bep*r_Der_Q,4)*r_Der_ss*(-3*P*r_Der_nn_2*(-1 + 8*Bep*r_Der_Q)*(3 + 880*pow(Bep,2)*pow(P,4) - 504*Bep*r_Der_Q + 7488*pow(Bep,2)*pow(r_Der_Q,2) - 29184*pow(Bep,3)*pow(r_Der_Q,3)) + nn*(3*r_Der_P_2*(-1 + 8*Bep*r_Der_Q)*(-7 - 2368*pow(Bep,2)*pow(P,4) + 184*Bep*r_Der_Q - 1600*pow(Bep,2)*pow(r_Der_Q,2) + 4608*pow(Bep,3)*pow(r_Der_Q,3)) + 16*Bep*P*(4*(-11 + 76*pow(Bep,2)*pow(P,4) + 264*Bep*r_Der_Q - 2112*pow(Bep,2)*pow(r_Der_Q,2) + 5632*pow(Bep,3)*pow(r_Der_Q,3))*r_Der_Q_3 + P*(3*P*(-17*Bepp*pow(P,2)*r_Der_Q - 141*Bep*pow(r_Der_Q,2) + 440*Bep*Bepp*pow(P,2)*pow(r_Der_Q,2) + 1128*pow(Bep,2)*pow(r_Der_Q,3)) + 8*(-2 - Bepp*pow(P,2) + 158*pow(Bep,2)*pow(P,4) + 4*Bep*(11 + 4*Bepp*pow(P,2))*r_Der_Q - 64*pow(Bep,2)*(5 + Bepp*pow(P,2))*pow(r_Der_Q,2) + 768*pow(Bep,3)*pow(r_Der_Q,3))*r_Der_ss_3)))))/12.;
      //
      // double p2denom = pow(1 - 8*Bep*r_Der_Q,8) + 72*Bep*P*pow(-1 + 8*Bep*r_Der_Q,7)*r_Der_ss + 2320*pow(Bep,2)*pow(P,2)*pow(1 - 8*Bep*r_Der_Q,6)*pow(r_Der_ss,2) + 43392*pow(Bep,3)*pow(P,3)*pow(-1 + 8*Bep*r_Der_Q,5)*pow(r_Der_ss,3) + 64*pow(Bep,2)*pow(1 - 8*Bep*r_Der_Q,4)*(-3 + 8008*pow(Bep,2)*pow(P,4) + 72*Bep*r_Der_Q - 576*pow(Bep,2)*pow(r_Der_Q,2) + 1536*pow(Bep,3)*pow(r_Der_Q,3))*pow(r_Der_ss,4) + 768*pow(Bep,3)*P*pow(-1 + 8*Bep*r_Der_Q,3)*(-13 + 5072*pow(Bep,2)*pow(P,4) + 312*Bep*r_Der_Q - 2496*pow(Bep,2)*pow(r_Der_Q,2) + 6656*pow(Bep,3)*pow(r_Der_Q,3))*pow(r_Der_ss,5) + 12288*pow(Bep,4)*pow(P,2)*pow(1 - 8*Bep*r_Der_Q,2)*(-18 + 1507*pow(Bep,2)*pow(P,4) + 432*Bep*r_Der_Q - 3456*pow(Bep,2)*pow(r_Der_Q,2) + 9216*pow(Bep,3)*pow(r_Der_Q,3))*pow(r_Der_ss,6) + 61440*pow(Bep,5)*pow(P,3)*(-1 + 8*Bep*r_Der_Q)*(-43 + 816*pow(Bep,2)*pow(P,4) + 1032*Bep*r_Der_Q - 8256*pow(Bep,2)*pow(r_Der_Q,2) + 22016*pow(Bep,3)*pow(r_Der_Q,3))*pow(r_Der_ss,7) + 3072*pow(Bep,4)*(3 - 5792*pow(Bep,2)*pow(P,4) + 19200*pow(Bep,4)*pow(P,8) + 48*Bep*(-3 + 2896*pow(Bep,2)*pow(P,4))*r_Der_Q - 192*pow(Bep,2)*(-15 + 5792*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,2) + 2048*pow(Bep,3)*(-15 + 1448*pow(Bep,2)*pow(P,4))*pow(r_Der_Q,3) + 184320*pow(Bep,4)*pow(r_Der_Q,4) - 589824*pow(Bep,5)*pow(r_Der_Q,5) + 786432*pow(Bep,6)*pow(r_Der_Q,6))*pow(r_Der_ss,8) + 294912*pow(Bep,5)*P*pow(1 - 8*Bep*r_Der_Q,2)*(-1 + 216*pow(Bep,2)*pow(P,4) + 24*Bep*r_Der_Q - 192*pow(Bep,2)*pow(r_Der_Q,2) + 512*pow(Bep,3)*pow(r_Der_Q,3))*pow(r_Der_ss,9) + 1179648*pow(Bep,6)*pow(P,2)*(-1 + 8*Bep*r_Der_Q)*(-3 + 80*pow(Bep,2)*pow(P,4) + 72*Bep*r_Der_Q - 576*pow(Bep,2)*pow(r_Der_Q,2) + 1536*pow(Bep,3)*pow(r_Der_Q,3))*pow(r_Der_ss,10) + 18874368*pow(Bep,7)*pow(P,3)*pow(-1 + 8*Bep*r_Der_Q,3)*pow(r_Der_ss,11) + 37748736*pow(Bep,8)*pow(P,4)*pow(1 - 8*Bep*r_Der_Q,2)*pow(r_Der_ss,12);
      //
      //
      // dpdt[i] += (p2num/p2denom)*(r[i]*r[i]);

      dqdt[i] = (P*r_Der_nn_2 + nn*(r_Der_P_2 + 2.*r_Der_Q*r_Der_ss))*r[i];

      dphidt[i] = n_v.v[i]*p_v.v[i] + r[i]*r[i]*((P*r_Der_nn_2)/2. + (nn*(r_Der_P_2 + 2.*r_Der_Q*r_Der_ss))/2.);
     }


  for(int i = grid.exc_i + 1; i<nx-1; i++){
    double r_Der_nn= Dx_ptpc_2nd(n_v.v[i+1], n_v.v[i-1], dr[i]);
    double r_Der_ss= Dx_ptpc_2nd(s_v.v[i+1], s_v.v[i-1], dr[i]);
    double r_Der_P= Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
    double r_Der_Q= Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
    double Bep = beta_p(ls,lexp,mu, phi_v.v[i]);
    double Bepp = beta_pp(ls,lexp,mu, phi_v.v[i]);

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
void Evolve_scalar_field::generate_rhs_excised(const Grid_data grid,
  const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v,
  double &dsdt,
  vector<double> &dpdt,
  vector<double> &dqdt,
  vector<double> &dphidt)
  {

  //Field ordering s_v, p_v,q_v,phi_v
  assert(grid.exc_i >0);
  vector<double> r = grid.r;
  int nx = grid.nx;
  vector<double> dr = grid.dr;
  double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;



  // {
  //   int i = grid.exc_i;
  //
  //   double Bep=  beta_p(ls,lexp,mu, phi_v.v[i]);
  //   double Bepp= beta_pp(ls,lexp,mu, phi_v.v[i]);
  //
  //   double r_Der_nn= Dx_ptp0_2nd(n_v.v[i+2], n_v.v[i+1], n_v.v[i], dr[i]);
  //   double r_Der_ss= Dx_ptp0_2nd(s_v.v[i+2], s_v.v[i+1], s_v.v[i], dr[i]);
  //   double r_Der_P= Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
  //   double r_Der_Q=Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);
  //
  //   dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q, Bep, Bepp);
  //   dqdt[i] = rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //   dphidt[i] = rhs_phi(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //   dsdt = rhs_s_free(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q, Bep, Bepp);
  //
  //
  // }
  //
  // for(int i = grid.exc_i + 1; i<nx-1; i++){
  //   double r_Der_nn= Dx_ptpc_2nd(n_v.v[i+1], n_v.v[i-1], dr[i]);
  //   double r_Der_ss= Dx_ptpc_2nd(s_v.v[i+1], s_v.v[i-1], dr[i]);
  //   double r_Der_P= Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
  //   double r_Der_Q= Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
  //   double Bep = beta_p(ls,lexp,mu, phi_v.v[i]);
  //   double Bepp = beta_pp(ls,lexp,mu, phi_v.v[i]);
  //
  //   dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q,Bep,Bepp);
  //   dqdt[i] = rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //   dphidt[i] = rhs_phi(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
  //
  //
  // }
  for(int i =grid.exc_i; i< grid.exc_i+1;i++){

    double Bep=  beta_p(ls,lexp,mu, phi_v.v[i]);
    double Bepp= beta_pp(ls,lexp,mu, phi_v.v[i]);

    double r_Der_nn= Dx_ptp0_2nd(n_v.v[i+2], n_v.v[i+1], n_v.v[i], dr[i]);
    double r_Der_ss= Dx_ptp0_2nd(s_v.v[i+2], s_v.v[i+1], s_v.v[i], dr[i]);
    double r_Der_P= Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
    double r_Der_Q=Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);

    dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q, Bep, Bepp);
    dqdt[i] = rhs_q(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
    dphidt[i] = rhs_phi(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q);
    dsdt = rhs_s_free(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q, Bep, Bepp);


  }

  for(int i = grid.exc_i+1; i<nx-1; i++){
    double r_Der_nn= Dx_ptpc_2nd(n_v.v[i+1], n_v.v[i-1], dr[i]);
    double r_Der_ss= Dx_ptpc_2nd(s_v.v[i+1], s_v.v[i-1], dr[i]);
    double r_Der_P= Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
    double r_Der_Q= Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);
    double Bep = beta_p(ls,lexp,mu, phi_v.v[i]);
    double Bepp = beta_pp(ls,lexp,mu, phi_v.v[i]);

    dpdt[i] = rhs_p(r[i], n_v.v[i], r_Der_nn, s_v.v[i], r_Der_ss, p_v.v[i], r_Der_P, q_v.v[i], r_Der_Q,Bep,Bepp);
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
//==============================================================================
void Evolve_scalar_field::evolve(const Grid_data grid, const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
  Solve_metric_fields solve_metric_fields;

  if(grid.exc_i>0){
    vector<double> r = grid.r;
    int nx = grid.nx;
    double dt = grid.dt;
    vector<double> dr = grid.dr;
    int exc_i = grid.exc_i;

    vector<double> kp1(nx,0), kq1(nx,0), kphi1(nx,0);
    vector<double> kp2(nx,0), kq2(nx,0), kphi2(nx,0);
    vector<double> kp3(nx,0), kq3(nx,0), kphi3(nx,0);
    vector<double> kp4(nx,0), kq4(nx,0), kphi4(nx,0);
    double sk1 = 0, sk2 = 0, sk3 = 0, sk4 = 0;
    vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);
    double dsdt;

    generate_rhs_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dsdt, dpdt, dqdt, dphidt);

    for(int i =exc_i; i<nx; i++){
      kp1[i] = dt*dpdt[i];
      kq1[i] = dt*dqdt[i];
      kphi1[i] = dt*dphidt[i];

    }
    sk1 = dt*dsdt;
    {
      Field n_1("n_k1", "even", grid);
      Field s_1("s_k1", "odd", grid);
      Field p_1("p_k1", "even", grid);
      Field q_1("q_k1", "odd", grid);
      Field phi_1("phi_k1", "even", grid);

      for(int i = exc_i; i<nx; i++){
        n_1.v[i] = n_v.v[i];
        s_1.v[i] = s_v.v[i] ;
        p_1.v[i] = p_v.v[i] + 0.5*kp1[i];
        q_1.v[i] = q_v.v[i] + 0.5*kq1[i];
        phi_1.v[i] = phi_v.v[i] + 0.5*kphi1[i];
      }
      s_1.v[exc_i] = s_v.v[exc_i] + 0.5*sk1;
      solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1,phi_1);
      generate_rhs_excised(grid, n_1, s_1, p_1, q_1, phi_1, dsdt, dpdt, dqdt, dphidt);

      for(int i =0; i<nx; i++){
        kp2[i] = dt*dpdt[i];
        kq2[i] = dt*dqdt[i];
        kphi2[i] = dt*dphidt[i];

      }
      sk2 = dt*dsdt;
    }

      {
        Field n_2("n_k2", "even", grid);
        Field s_2("s_k2", "odd", grid);
        Field p_2("p_k2", "even", grid);
        Field q_2("q_k2", "odd", grid);
        Field phi_2("phi_k2", "even", grid);

        for(int i = exc_i; i<nx; i++){
          n_2.v[i] = n_v.v[i];
          s_2.v[i] = s_v.v[i] ;
          p_2.v[i] = p_v.v[i] + 0.5*kp2[i];
          q_2.v[i] = q_v.v[i] + 0.5*kq2[i];
          phi_2.v[i] = phi_v.v[i] + 0.5*kphi2[i];
        }
        s_2.v[exc_i] = s_v.v[exc_i] + 0.5*sk2;
        solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2,phi_2);
        generate_rhs_excised(grid, n_2, s_2, p_2, q_2, phi_2, dsdt, dpdt, dqdt, dphidt);

        for(int i =exc_i; i<nx; i++){
          kp3[i] = dt*dpdt[i];
          kq3[i] = dt*dqdt[i];
          kphi3[i] = dt*dphidt[i];

        }
        sk3 = dt*dsdt;



      }
      {
        Field n_3("n_k3", "even", grid);
        Field s_3("s_k3", "odd", grid);
        Field p_3("p_k3", "even", grid);
        Field q_3("q_k3", "odd", grid);
        Field phi_3("phi_k3", "even", grid);

        for(int i = exc_i; i<nx; i++){
          n_3.v[i] = n_v.v[i];
          s_3.v[i] = s_v.v[i];
          p_3.v[i] = p_v.v[i] + kp3[i];
          q_3.v[i] = q_v.v[i] + kq3[i];
          phi_3.v[i] = phi_v.v[i] + kphi3[i];
        }
        s_3.v[exc_i] = s_v.v[exc_i] + sk3;
        solve_metric_fields.solve(grid, n_3, s_3, p_3, q_3,phi_3);
        generate_rhs_excised(grid, n_3, s_3, p_3, q_3, phi_3, dsdt, dpdt, dqdt, dphidt);

        for(int i =exc_i; i<nx; i++){
          kp4[i] = dt*dpdt[i];
          kq4[i] = dt*dqdt[i];
          kphi4[i] = dt*dphidt[i];

        }
        sk4 = dt*dsdt;



      }
      for(int i=exc_i; i<nx; i++){
        p_v.v[i] += (1./6.)*kp1[i] + (1./3.)*kp2[i] + (1./3.)*kp3[i] + (1./6.)*kp4[i];
        q_v.v[i] += (1./6.)*kq1[i] + (1./3.)*kq2[i] + (1./3.)*kq3[i] + (1./6.)*kq4[i];
        phi_v.v[i] += (1./6.)*kphi1[i] + (1./3.)*kphi2[i] + (1./3.)*kphi3[i] + (1./6.)*kphi4[i];
      }
      for(int i=0; i<exc_i; i++){
        p_v.v[i] = 0.;
        q_v.v[i] = 0.;
        phi_v.v[i] = 0.;
        s_v.v[i] = 0.;
      }
      s_v.v[exc_i] += (1./6.)*sk1 + (1./3.)*sk2 + (1./3.)*sk3 + (1./6.)*sk4;

      p_v.check_isfinite(grid.t_evolve);
      q_v.check_isfinite(grid.t_evolve);
      phi_v.check_isfinite(grid.t_evolve);

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

    for(int i =0; i<nx; i++){
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

      for(int i = 0; i<nx; i++){
        n_1.v[i] = n_v.v[i];
        s_1.v[i] = s_v.v[i];
        p_1.v[i] = p_v.v[i] + 0.5*kp1[i];
        q_1.v[i] = q_v.v[i] + 0.5*kq1[i];
        phi_1.v[i] = phi_v.v[i] + 0.5*kphi1[i];
      }
      solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1, phi_1);
      generate_rhs_non_excised(grid, n_1, s_1, p_1, q_1, phi_1, dpdt, dqdt, dphidt);

      for(int i =0; i<nx; i++){
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

        for(int i = 0; i<nx; i++){
          n_2.v[i] = n_v.v[i];
          s_2.v[i] = s_v.v[i];
          p_2.v[i] = p_v.v[i] + 0.5*kp2[i];
          q_2.v[i] = q_v.v[i] + 0.5*kq2[i];
          phi_2.v[i] = phi_v.v[i] + 0.5*kphi2[i];
        }
        solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2, phi_2);
        generate_rhs_non_excised(grid, n_2, s_2, p_2, q_2, phi_2, dpdt, dqdt, dphidt);

        for(int i =0; i<nx; i++){
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

        for(int i = 0; i<nx; i++){
          n_3.v[i] = n_v.v[i];
          s_3.v[i] = s_v.v[i];
          p_3.v[i] = p_v.v[i] + kp3[i];
          q_3.v[i] = q_v.v[i] + kq3[i];
          phi_3.v[i] = phi_v.v[i] + kphi3[i];
        }
        solve_metric_fields.solve(grid, n_3, s_3, p_3, q_3, phi_3);
        generate_rhs_non_excised(grid, n_3, s_3, p_3, q_3, phi_3, dpdt, dqdt, dphidt);

        for(int i =0; i<nx; i++){
          kp4[i] = dt*dpdt[i];
          kq4[i] = dt*dqdt[i];
          kphi4[i] = dt*dphidt[i];

        }



      }
      for(int i=0; i<nx; i++){
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
//==============================================================================
//RK3 SSRPK2
// void Evolve_scalar_field::evolve(const Grid_data grid, const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
//   Solve_metric_fields solve_metric_fields;
//
//   if(grid.exc_i>0){
//     vector<double> r = grid.r;
//     int nx = grid.nx;
//     double dt = grid.dt;
//     vector<double> dr = grid.dr;
//     int exc_i = grid.exc_i;
//
//     vector<double> kp1(nx,0), kq1(nx,0), kphi1(nx,0);
//     vector<double> kp2(nx,0), kq2(nx,0), kphi2(nx,0);
//     vector<double> kp3(nx,0), kq3(nx,0), kphi3(nx,0);
//     double sk1 = 0, sk2 = 0, sk3 = 0;
//     vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);
//     double dsdt;
//
//     generate_rhs_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dsdt, dpdt, dqdt, dphidt);
//
//     for(int i =exc_i; i<nx; i++){
//       kp1[i] = dt*dpdt[i];
//       kq1[i] = dt*dqdt[i];
//       kphi1[i] = dt*dphidt[i];
//
//     }
//     sk1 = dt*dsdt;
//     {
//       Field n_1("n_k1", "even", grid);
//       Field s_1("s_k1", "odd", grid);
//       Field p_1("p_k1", "even", grid);
//       Field q_1("q_k1", "odd", grid);
//       Field phi_1("phi_k1", "even", grid);
//
//       for(int i = exc_i; i<nx; i++){
//         n_1.v[i] = n_v.v[i];
//         s_1.v[i] = s_v.v[i] ;
//         p_1.v[i] = p_v.v[i] + kp1[i];
//         q_1.v[i] = q_v.v[i] + kq1[i];
//         phi_1.v[i] = phi_v.v[i] + kphi1[i];
//       }
//       s_1.v[exc_i] = s_v.v[exc_i] + sk1;
//       solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1,phi_1);
//       generate_rhs_excised(grid, n_1, s_1, p_1, q_1, phi_1, dsdt, dpdt, dqdt, dphidt);
//
//       for(int i =0; i<nx-1; i++){
//         kp2[i] = dt*dpdt[i];
//         kq2[i] = dt*dqdt[i];
//         kphi2[i] = dt*dphidt[i];
//
//       }
//       sk2 = dt*dsdt;
//     }
//
//       {
//         Field n_2("n_k2", "even", grid);
//         Field s_2("s_k2", "odd", grid);
//         Field p_2("p_k2", "even", grid);
//         Field q_2("q_k2", "odd", grid);
//         Field phi_2("phi_k2", "even", grid);
//
//         for(int i = exc_i; i<nx; i++){
//           n_2.v[i] = n_v.v[i];
//           s_2.v[i] = s_v.v[i] ;
//           p_2.v[i] = p_v.v[i] + 0.25*kp1[i] + 0.25*kp2[i];
//           q_2.v[i] = q_v.v[i] + 0.25*kq1[i] + 0.25*kq2[i];
//           phi_2.v[i] = phi_v.v[i] + 0.25*kphi1[i] + 0.25*kphi2[i];
//         }
//         s_2.v[exc_i] = s_v.v[exc_i] + 0.25*sk1 + 0.25*sk2;
//         solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2,phi_2);
//         generate_rhs_excised(grid, n_2, s_2, p_2, q_2, phi_2, dsdt, dpdt, dqdt, dphidt);
//
//         for(int i =exc_i; i<nx-1; i++){
//           kp3[i] = dt*dpdt[i];
//           kq3[i] = dt*dqdt[i];
//           kphi3[i] = dt*dphidt[i];
//
//         }
//         sk3 = dt*dsdt;
//
//
//
//       }
//
//       for(int i=exc_i; i<nx; i++){
//         p_v.v[i] += (1./6.)*kp1[i] + (1./6.)*kp2[i] + (2./3.)*kp3[i] ;
//         q_v.v[i] += (1./6.)*kq1[i] + (1./6.)*kq2[i] + (2./3.)*kq3[i] ;
//         phi_v.v[i] += (1./6.)*kphi1[i] + (1./6.)*kphi2[i] + (2./3.)*kphi3[i] ;
//       }
//       for(int i=0; i<exc_i; i++){
//         p_v.v[i] = p_v.v[exc_i];
//         q_v.v[i] = q_v.v[exc_i];
//         phi_v.v[i] = phi_v.v[exc_i];
//       }
//       s_v.v[exc_i] += (1./6.)*sk1 + (1./6.)*sk2 + (2./3.)*sk3 ;
//
//       p_v.check_isfinite(grid.t_evolve);
//       q_v.check_isfinite(grid.t_evolve);
//       phi_v.check_isfinite(grid.t_evolve);
//
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
//     vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);
//
//     generate_rhs_non_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dpdt, dqdt, dphidt);
//
//     for(int i =0; i<nx; i++){
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
//       for(int i = 0; i<nx; i++){
//         n_1.v[i] = n_v.v[i];
//         s_1.v[i] = s_v.v[i];
//         p_1.v[i] = p_v.v[i] + kp1[i];
//         q_1.v[i] = q_v.v[i] + kq1[i];
//         phi_1.v[i] = phi_v.v[i] + kphi1[i];
//       }
//       solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1, phi_1);
//       generate_rhs_non_excised(grid, n_1, s_1, p_1, q_1, phi_1, dpdt, dqdt, dphidt);
//
//       for(int i =0; i<nx; i++){
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
//         for(int i = 0; i<nx; i++){
//           n_2.v[i] = n_v.v[i];
//           s_2.v[i] = s_v.v[i];
//           p_2.v[i] = p_v.v[i] + 0.25*kp1[i] + 0.25*kp2[i];
//           q_2.v[i] = q_v.v[i] + 0.25*kq1[i] + 0.25*kq2[i];
//           phi_2.v[i] = phi_v.v[i] + 0.25*kphi1[i] + 0.25*kphi2[i];
//         }
//         solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2, phi_2);
//         generate_rhs_non_excised(grid, n_2, s_2, p_2, q_2, phi_2, dpdt, dqdt, dphidt);
//
//         for(int i =0; i<nx; i++){
//           kp3[i] = dt*dpdt[i];
//           kq3[i] = dt*dqdt[i];
//           kphi3[i] = dt*dphidt[i];
//
//         }
//
//
//
//       }
//       for(int i=0; i<nx; i++){
//         p_v.v[i] += (1./6.)*kp1[i] + (1./6.)*kp2[i] + (2./3.)*kp3[i] ;
//         q_v.v[i] += (1./6.)*kq1[i] + (1./6.)*kq2[i] + (2./3.)*kq3[i] ;
//         phi_v.v[i] += (1./6.)*kphi1[i] + (1./6.)*kphi2[i] + (2./3.)*kphi3[i] ;
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
//==========================================================================================================================
//==============================================================================
// void Evolve_scalar_field::evolve(const Grid_data grid, const Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
//   Solve_metric_fields solve_metric_fields;
//
//   if(grid.exc_i>0){
//     vector<double> r = grid.r;
//     int nx = grid.nx;
//     double dt = grid.dt;
//     vector<double> dr = grid.dr;
//     int exc_i = grid.exc_i;
//
//     vector<double> kp1(nx,0), kq1(nx,0), kphi1(nx,0);
//     vector<double> kp2(nx,0), kq2(nx,0), kphi2(nx,0);
//     vector<double> kp3(nx,0), kq3(nx,0), kphi3(nx,0);
//     vector<double> kp4(nx,0), kq4(nx,0), kphi4(nx,0);
//     double sk1 = 0, sk2 = 0, sk3 = 0, sk4 = 0;
//     vector<double> dpdt(nx,0), dqdt(nx,0), dphidt(nx,0);
//     double dsdt;
//
//     generate_rhs_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dsdt, dpdt, dqdt, dphidt);
//
//     for(int i =exc_i; i<nx-1; i++){
//       kp1[i] = dt*dpdt[i];
//       kq1[i] = dt*dqdt[i];
//       kphi1[i] = dt*dphidt[i];
//
//     }
//     sk1 = dt*dsdt;
//     {
//       Field n_1("n_k1", "even", grid);
//       Field s_1("s_k1", "odd", grid);
//       Field p_1("p_k1", "even", grid);
//       Field q_1("q_k1", "odd", grid);
//       Field phi_1("phi_k1", "even", grid);
//
//       for(int i = exc_i; i<nx-1; i++){
//         n_1.v[i] = n_v.v[i];
//         s_1.v[i] = s_v.v[i] ;
//         p_1.v[i] = p_v.v[i] + 0.5*kp1[i];
//         q_1.v[i] = q_v.v[i] + 0.5*kq1[i];
//         phi_1.v[i] = phi_v.v[i] + 0.5*kphi1[i];
//       }
//       s_1.v[exc_i] = s_v.v[exc_i] + 0.5*sk1;
//       solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1,phi_1);
//       generate_rhs_excised(grid, n_1, s_1, p_1, q_1, phi_1, dsdt, dpdt, dqdt, dphidt);
//
//       for(int i =0; i<nx-1; i++){
//         kp2[i] = dt*dpdt[i];
//         kq2[i] = dt*dqdt[i];
//         kphi2[i] = dt*dphidt[i];
//
//       }
//       sk2 = dt*dsdt;
//     }
//
//       {
//         Field n_2("n_k2", "even", grid);
//         Field s_2("s_k2", "odd", grid);
//         Field p_2("p_k2", "even", grid);
//         Field q_2("q_k2", "odd", grid);
//         Field phi_2("phi_k2", "even", grid);
//
//         for(int i = exc_i; i<nx-1; i++){
//           n_2.v[i] = n_v.v[i];
//           s_2.v[i] = s_v.v[i] ;
//           p_2.v[i] = p_v.v[i] + 0.5*kp2[i];
//           q_2.v[i] = q_v.v[i] + 0.5*kq2[i];
//           phi_2.v[i] = phi_v.v[i] + 0.5*kphi2[i];
//         }
//         s_2.v[exc_i] = s_v.v[exc_i] + 0.5*sk2;
//         solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2,phi_2);
//         generate_rhs_excised(grid, n_2, s_2, p_2, q_2, phi_2, dsdt, dpdt, dqdt, dphidt);
//
//         for(int i =exc_i; i<nx-1; i++){
//           kp3[i] = dt*dpdt[i];
//           kq3[i] = dt*dqdt[i];
//           kphi3[i] = dt*dphidt[i];
//
//         }
//         sk3 = dt*dsdt;
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
//         for(int i = exc_i; i<nx-1; i++){
//           n_3.v[i] = n_v.v[i];
//           s_3.v[i] = s_v.v[i];
//           p_3.v[i] = p_v.v[i] + kp3[i];
//           q_3.v[i] = q_v.v[i] + kq3[i];
//           phi_3.v[i] = phi_v.v[i] + kphi3[i];
//         }
//         s_3.v[exc_i] = s_v.v[exc_i] + sk3;
//         solve_metric_fields.solve(grid, n_3, s_3, p_3, q_3,phi_3);
//         generate_rhs_excised(grid, n_3, s_3, p_3, q_3, phi_3, dsdt, dpdt, dqdt, dphidt);
//
//         for(int i =exc_i; i<nx-1; i++){
//           kp4[i] = dt*dpdt[i];
//           kq4[i] = dt*dqdt[i];
//           kphi4[i] = dt*dphidt[i];
//
//         }
//         sk4 = dt*dsdt;
//
//
//
//       }
//       for(int i=exc_i; i<nx-1; i++){
//         p_v.v[i] += (1./6.)*kp1[i] + (1./3.)*kp2[i] + (1./3.)*kp3[i] + (1./6.)*kp4[i];
//         q_v.v[i] += (1./6.)*kq1[i] + (1./3.)*kq2[i] + (1./3.)*kq3[i] + (1./6.)*kq4[i];
//         phi_v.v[i] += (1./6.)*kphi1[i] + (1./3.)*kphi2[i] + (1./3.)*kphi3[i] + (1./6.)*kphi4[i];
//       }
//       for(int i=0; i<exc_i; i++){
//         p_v.v[i] = p_v.v[exc_i];
//         q_v.v[i] = q_v.v[exc_i];
//         phi_v.v[i] = phi_v.v[exc_i];
//       }
//       s_v.v[exc_i] += (1./6.)*sk1 + (1./3.)*sk2 + (1./3.)*sk3 + (1./6.)*sk4;
//
//       p_v.check_isfinite(grid.t_evolve);
//       q_v.check_isfinite(grid.t_evolve);
//       phi_v.check_isfinite(grid.t_evolve);
//
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
//     generate_rhs_non_excised(grid ,n_v, s_v, p_v, q_v, phi_v, dpdt, dqdt, dphidt);
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
//       solve_metric_fields.solve(grid, n_1, s_1, p_1, q_1, phi_1);
//       generate_rhs_non_excised(grid, n_1, s_1, p_1, q_1, phi_1, dpdt, dqdt, dphidt);
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
//         solve_metric_fields.solve(grid, n_2, s_2, p_2, q_2, phi_2);
//         generate_rhs_non_excised(grid, n_2, s_2, p_2, q_2, phi_2, dpdt, dqdt, dphidt);
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
//         solve_metric_fields.solve(grid, n_3, s_3, p_3, q_3, phi_3);
//         generate_rhs_non_excised(grid, n_3, s_3, p_3, q_3, phi_3, dpdt, dqdt, dphidt);
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
//         p_v.v[i] += (1./6.)*kp1[i] + (1./3.)*kp2[i] + (1./3.)*kp3[i] + (1./6.)*kp4[i];
//         q_v.v[i] += (1./6.)*kq1[i] + (1./3.)*kq2[i] + (1./3.)*kq3[i] + (1./6.)*kq4[i];
//         phi_v.v[i] += (1./6.)*kphi1[i] + (1./3.)*kphi2[i] + (1./3.)*kphi3[i] + (1./6.)*kphi4[i];
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
