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
    if((fabs(ss)<1e-8)&& r<12.0){
      return 2*Qr*nn + r_Der_Q*nn + r*r_Der_P*ssr*nn + r_Der_ss*nn*P + 2*ssr*nn*P + r_Der_nn*(Qr*r + r*ssr*P);
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
   double eps= grid.dissipation;
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
      int i = grid.exc_i;

      double nn = n_v.v[i];
      double r_Der_nn = 0.;

      double ss = 0.;
      double r_Der_ss = Dx_ptpc_2nd(s_v.v[i+1], -s_v.v[i+1], dr[i]);

      double P = p_v.v[i];
      double r_Der_P = 0.;

      double Q = 0.;
      double r_Der_Q = Dx_ptpc_2nd(q_v.v[i+1], -q_v.v[i+1], dr[i]);

      double Bep = beta_p(ls,lexp,mu, phi_v.v[i]);
      double Bepp = beta_pp(ls,lexp,mu, phi_v.v[i]);

      double rr_Der_nn = Dx_2_ptpc_2nd(n_v.v[i+1], n_v.v[i], n_v.v[i+1], dr[i]);
      double rr_Der_P =Dx_2_ptpc_2nd(p_v.v[i+1], p_v.v[i], p_v.v[i+1], dr[i]);

      dpdt[i] = nn*(-pow(P,3) + 16*Bep*pow(P,3)*r_Der_Q - 64*pow(Bep,2)*pow(P,3)*pow(r_Der_Q,2) + 20*Bep*pow(P,4)*r_Der_ss - 12*r_Der_Q*r_Der_ss - 160*pow(Bep,2)*pow(P,4)*r_Der_Q*r_Der_ss + 288*Bep*pow(r_Der_Q,2)*r_Der_ss - 2304*pow(Bep,2)*pow(r_Der_Q,3)*r_Der_ss + 6144*pow(Bep,3)*pow(r_Der_Q,4)*r_Der_ss - 6*P*pow(r_Der_ss,2) - 80*pow(Bep,2)*pow(P,5)*pow(r_Der_ss,2) + 480*Bep*P*r_Der_Q*pow(r_Der_ss,2) - 6528*pow(Bep,2)*P*pow(r_Der_Q,2)*pow(r_Der_ss,2) + 24576*pow(Bep,3)*P*pow(r_Der_Q,3)*pow(r_Der_ss,2) + 208*Bep*pow(P,2)*pow(r_Der_ss,3) + 128*Bep*Bepp*pow(P,4)*pow(r_Der_ss,3) - 6272*pow(Bep,2)*pow(P,2)*r_Der_Q*pow(r_Der_ss,3) - 1024*pow(Bep,2)*Bepp*pow(P,4)*r_Der_Q*pow(r_Der_ss,3) + 36864*pow(Bep,3)*pow(P,2)*pow(r_Der_Q,2)*pow(r_Der_ss,3) - 2208*pow(Bep,2)*pow(P,3)*pow(r_Der_ss,4) - 768*pow(Bep,2)*Bepp*pow(P,5)*pow(r_Der_ss,4) + 26112*pow(Bep,3)*pow(P,3)*r_Der_Q*pow(r_Der_ss,4) - 48*Bep*pow(r_Der_ss,5) - 384*Bep*Bepp*pow(P,2)*pow(r_Der_ss,5) + 7680*pow(Bep,3)*pow(P,4)*pow(r_Der_ss,5) + 768*pow(Bep,2)*r_Der_Q*pow(r_Der_ss,5) + 6144*pow(Bep,2)*Bepp*pow(P,2)*r_Der_Q*pow(r_Der_ss,5) - 3072*pow(Bep,3)*pow(r_Der_Q,2)*pow(r_Der_ss,5) - 24576*pow(Bep,3)*Bepp*pow(P,2)*pow(r_Der_Q,2)*pow(r_Der_ss,5) + 384*pow(Bep,2)*P*pow(r_Der_ss,6) + 3072*pow(Bep,2)*Bepp*pow(P,3)*pow(r_Der_ss,6) - 3072*pow(Bep,3)*P*r_Der_Q*pow(r_Der_ss,6) - 24576*pow(Bep,3)*Bepp*pow(P,3)*r_Der_Q*pow(r_Der_ss,6))
      ;

      dpdt[i] /= 4*r_Der_ss*(-1 + 24*Bep*r_Der_Q - 192*pow(Bep,2)*pow(r_Der_Q,2) + 512*pow(Bep,3)*pow(r_Der_Q,3) + 28*Bep*P*r_Der_ss - 448*pow(Bep,2)*P*r_Der_Q*r_Der_ss + 1792*pow(Bep,3)*P*pow(r_Der_Q,2)*r_Der_ss - 288*pow(Bep,2)*pow(P,2)*pow(r_Der_ss,2) + 2304*pow(Bep,3)*pow(P,2)*r_Der_Q*pow(r_Der_ss,2) + 960*pow(Bep,3)*pow(P,3)*pow(r_Der_ss,3) + 96*pow(Bep,2)*pow(r_Der_ss,4) - 1536*pow(Bep,3)*r_Der_Q*pow(r_Der_ss,4) + 6144*pow(Bep,4)*pow(r_Der_Q,2)*pow(r_Der_ss,4) - 768*pow(Bep,3)*P*pow(r_Der_ss,5) + 6144*pow(Bep,4)*P*r_Der_Q*pow(r_Der_ss,5))
      ;

      dqdt[i] = 0.;

      dphidt[i] = nn*P;

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
  if(grid.bh_start==0){
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
  else{

    for(int i =grid.exc_i; i< grid.ah_index-3;i++){

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

    for(int i = grid.ah_index-3; i<nx-1; i++){
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
