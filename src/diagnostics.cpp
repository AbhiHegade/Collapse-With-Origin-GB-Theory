#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
using std::pow;
using std::fabs;
using std::exp;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include<cassert>
#include <algorithm>
using std::min_element;

#include "field.hpp"
#include "grid_data.hpp"
#include "diagnostics.hpp"
#include "compute_potentials.hpp"
#include "fd_stencils.hpp"
//==============================================================================
Diagnostics::Diagnostics(){

}
//==============================================================================
Diagnostics::~Diagnostics(void){

}
//==============================================================================
void Diagnostics::find_abs_min(const vector<double> &v,
  double &min_elem,
  int &index,
  const double ref_val,
  const int start_index )
  {
    int len = v.size();

    vector<double> v_abs(len,1);

    for(int i=start_index; i<len;i++){
      v_abs[i] = fabs(ref_val - v[i]);
    }

    auto minret = min_element(v_abs.begin(), v_abs.end());

    min_elem = *minret;

    index = minret - v_abs.begin();

    return;
}
//==============================================================================
void Diagnostics::find_outer_most_index(const vector<double> &v, int &elem, const int start_index)
{

  double ref_val = 1.0;
  double tol = 1e-2;
  int index_val = -1;
  int len = 3*(v.size())/4;
  for(int i = start_index; i<len; i++){
    if((fabs(ref_val - v[i])<=tol))
    {
      index_val = i;
      // cout<<v[i]<<"\t"<<index_val<<endl;
    }
  }
  elem = index_val;
  return;

}
//==============================================================================
void Diagnostics::find_apparent_horizon(Grid_data &grid, Field &s_v)
{
  int index = 0;
  find_outer_most_index(s_v.v,index,1);
  // cout<<"index = "<<index<<endl;
  if(index>=0){
    if(grid.exc_i>0){
      int ex_int = grid.ex_ratio*index;
      int indexby2 = ( ex_int ==0 ) ? (index/2) : (ex_int);
      int new_exc_i = (indexby2>grid.exc_i) ? indexby2 : grid.exc_i;

      if (index==0){
          cout<<"Apparent horizon at the origin."<<endl;
          std::exit(0);
        }
      else if (new_exc_i == grid.exc_i){
          grid.ah_index = index;
        }
      else{
          grid.ah_index = index;
          cout<<"Found apparent horizon at i = "<<index<<" , "<<"r = "<<grid.r[index]<<" , "<< "t = "<< grid.t_evolve<<endl;
          cout<<"Previous excision point at i = "<<grid.exc_i<<" , "<<"r = "<<grid.r[grid.exc_i]<<endl;
          cout<<"Updating excision point to i = "<<new_exc_i<<", r = "<<grid.r[new_exc_i]<<endl;
          grid.exc_i = new_exc_i;
          cout<<"done."<<endl;
        }

    }

    else{
      assert(grid.exc_i==0);
      assert(grid.ah_index ==0);
      int ex_int = grid.ex_ratio*index;
      int indexby2 = ( ex_int ==0 ) ? (index/2) : (ex_int);
      if (index==0){
          cout<<"Apparent Horizon at the origin."<<endl;
          std::exit(0);
        }
      else{
          // int indexby2 = (index%2 ==0 ) ? (index/2) : (index+1)/2;
          grid.ah_index = index;
          // int indexby2 = (index/pt ==0) ? (index/(pt-1)) : (index/pt);

          cout<<"Found apparent horizon at i = "<<index<<" , "<<"r = "<<grid.r[index]<<" , "<< "t = "<< grid.t_evolve<<endl;
          cout<<"Previous excision point at i = "<<grid.exc_i<<" , "<<"r = "<<grid.r[grid.exc_i]<<endl;
          cout<<"Updating excision point to i = "<<indexby2<<", r = "<<grid.r[indexby2]<<", ";
          grid.exc_i = indexby2;
          cout<<"done."<<endl;
        }
      }
    }
}
//==============================================================================
int Diagnostics::compute_radial_characteristic(double r,
double nn, double r_Der_nn,
double ss, double r_Der_ss,
double P, double r_Der_P,
double Q, double r_Der_Q,
double Bep, double Bepp,
double &ingoing_c, double &outgoing_c)
{

  double Qr= Q/r;
  double ssr= ss/r;

  // dP/dt .
  double t_Der_P= ((-64*pow(r_Der_nn,2)*(r - 8*Bep*Qr*r - 8*Bep*r*ssr*P)*(pow(Bep,2)*Qr*pow(r,5)*pow(ssr,4) + pow(Bep,2)*pow(r,3)*pow(ssr,3)*P))/nn - 64*pow(r_Der_ss,2)*(pow(Bep,2)*Qr*pow(r,4)*pow(ssr,2)*nn - 8*pow(Bep,3)*pow(Qr,2)*pow(r,4)*pow(ssr,2)*nn - 4*pow(Bep,3)*Qr*pow(r,4)*pow(ssr,3)*nn*P) - r_Der_P*(-(pow(r,5)*ssr*nn) + 16*Bep*Qr*pow(r,5)*ssr*nn - 64*pow(Bep,2)*pow(Qr,2)*pow(r,5)*ssr*nn - 32*pow(Bep,2)*pow(r,3)*pow(ssr,3)*nn + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,5)*pow(ssr,3)*nn + 32*pow(Bep,2)*pow(r,5)*pow(ssr,5)*nn - 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,7)*pow(ssr,5)*nn - 256*r_Der_Q*(-(pow(Bep,3)*pow(r,3)*pow(ssr,3)*nn) + pow(Bep,3)*pow(r,5)*pow(ssr,5)*nn) + 16*Bep*pow(r,5)*pow(ssr,2)*nn*P - 128*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,2)*nn*P - 64*pow(Bep,2)*pow(r,5)*pow(ssr,3)*nn*pow(P,2)) - r_Der_Q*(-(pow(r,4)*nn) + 16*Bep*Qr*pow(r,4)*nn - 64*pow(Bep,2)*pow(Qr,2)*pow(r,4)*nn + 16*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn + 32*pow(Bep,2)*pow(r,4)*pow(ssr,4)*nn + 16*Bep*pow(r,4)*ssr*nn*P - 96*pow(Bep,2)*Qr*pow(r,4)*ssr*nn*P + 256*pow(Bep,2)*Bepp*Qr*pow(r,4)*pow(ssr,3)*nn*P - 48*pow(Bep,2)*pow(r,4)*pow(ssr,2)*nn*pow(P,2) + 256*pow(Bep,2)*Bepp*pow(r,4)*pow(ssr,4)*nn*pow(P,2)) + (8*Qr*pow(r,4)*nn - 128*Bep*pow(Qr,2)*pow(r,4)*nn + 512*pow(Bep,2)*pow(Qr,3)*pow(r,4)*nn + 24*Bep*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn - 128*pow(Bep,2)*pow(Qr,3)*pow(r,6)*pow(ssr,2)*nn - 64*Bep*Bepp*pow(Qr,4)*pow(r,8)*pow(ssr,2)*nn + 16*Bep*pow(r,4)*pow(ssr,4)*nn - 128*Bep*Bepp*pow(Qr,2)*pow(r,6)*pow(ssr,4)*nn + 8*pow(r,4)*ssr*nn*P - 240*Bep*Qr*pow(r,4)*ssr*nn*P + 1536*pow(Bep,2)*pow(Qr,2)*pow(r,4)*ssr*nn*P - 128*Bep*Bepp*pow(Qr,3)*pow(r,6)*ssr*nn*P + 128*Bep*Bepp*Qr*pow(r,4)*pow(ssr,3)*nn*P - 128*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(ssr,3)*nn*P - 1024*Bep*pow(Bepp,2)*pow(Qr,3)*pow(r,6)*pow(ssr,3)*nn*P - 136*Bep*pow(r,4)*pow(ssr,2)*nn*pow(P,2) + 1664*pow(Bep,2)*Qr*pow(r,4)*pow(ssr,2)*nn*pow(P,2) - 64*Bep*Bepp*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn*pow(P,2) + 128*Bep*Bepp*pow(r,4)*pow(ssr,4)*nn*pow(P,2) - 1024*Bep*pow(Bepp,2)*pow(Qr,2)*pow(r,6)*pow(ssr,4)*nn*pow(P,2) + 640*pow(Bep,2)*pow(r,4)*pow(ssr,3)*nn*pow(P,3))/4. - r_Der_ss*(4*Bep*pow(Qr,2)*pow(r,6)*ssr*nn - 32*pow(Bep,2)*pow(Qr,3)*pow(r,6)*ssr*nn + 16*Bep*pow(r,4)*pow(ssr,3)*nn - 160*pow(Bep,2)*Qr*pow(r,4)*pow(ssr,3)*nn + 256*pow(Bep,2)*Bepp*pow(Qr,3)*pow(r,6)*pow(ssr,3)*nn + 256*pow(Bep,3)*Qr*pow(r,4)*r_Der_Q*pow(ssr,3)*nn - pow(r,4)*nn*P + 24*Bep*Qr*pow(r,4)*nn*P - 128*pow(Bep,2)*pow(Qr,2)*pow(r,4)*nn*P + 64*Bep*Bepp*Qr*pow(r,4)*pow(ssr,2)*nn*P - 512*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,4)*pow(ssr,2)*nn*P - 16*pow(Bep,2)*pow(Qr,2)*pow(r,6)*pow(ssr,2)*nn*P - 96*pow(Bep,2)*pow(r,4)*pow(ssr,4)*nn*P + 20*Bep*pow(r,4)*ssr*nn*pow(P,2) - 192*pow(Bep,2)*Qr*pow(r,4)*ssr*nn*pow(P,2) + 128*Bep*Bepp*pow(r,4)*pow(ssr,3)*nn*pow(P,2) - 1280*pow(Bep,2)*Bepp*Qr*pow(r,4)*pow(ssr,3)*nn*pow(P,2) - 80*pow(Bep,2)*pow(r,4)*pow(ssr,2)*nn*pow(P,3) - 768*pow(Bep,2)*Bepp*pow(r,4)*pow(ssr,4)*nn*pow(P,3) - 64*r_Der_P*(-(pow(Bep,2)*pow(r,3)*pow(ssr,2)*nn) + 8*pow(Bep,3)*Qr*pow(r,3)*pow(ssr,2)*nn + 2*pow(Bep,2)*pow(r,5)*pow(ssr,4)*nn - 16*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4)*nn + 4*pow(Bep,3)*pow(r,3)*pow(ssr,3)*nn*P - 12*pow(Bep,3)*pow(r,5)*pow(ssr,5)*nn*P)) - r_Der_nn*(-(Qr*pow(r,5)) + 16*Bep*pow(Qr,2)*pow(r,5) - 64*pow(Bep,2)*pow(Qr,3)*pow(r,5) - 8*Bep*pow(r,3)*pow(ssr,2) + 64*pow(Bep,2)*Qr*pow(r,3)*pow(ssr,2) + 64*Bep*Bepp*pow(Qr,2)*pow(r,5)*pow(ssr,2) - 512*pow(Bep,2)*Bepp*pow(Qr,3)*pow(r,5)*pow(ssr,2) + 4*Bep*pow(Qr,2)*pow(r,7)*pow(ssr,2) - 32*pow(Bep,2)*pow(Qr,3)*pow(r,7)*pow(ssr,2) + 16*Bep*pow(r,5)*pow(ssr,4) - 128*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,4) - pow(r,5)*ssr*P + 40*Bep*Qr*pow(r,5)*ssr*P - 256*pow(Bep,2)*pow(Qr,2)*pow(r,5)*ssr*P + 64*pow(Bep,2)*pow(r,3)*pow(ssr,3)*P + 192*Bep*Bepp*Qr*pow(r,5)*pow(ssr,3)*P - 2048*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,5)*pow(ssr,3)*P - 16*pow(Bep,2)*pow(Qr,2)*pow(r,7)*pow(ssr,3)*P - 96*pow(Bep,2)*pow(r,5)*pow(ssr,5)*P + 20*Bep*pow(r,5)*pow(ssr,2)*pow(P,2) - 256*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,2)*pow(P,2) + 128*Bep*Bepp*pow(r,5)*pow(ssr,4)*pow(P,2) - 2304*pow(Bep,2)*Bepp*Qr*pow(r,5)*pow(ssr,4)*pow(P,2) - 80*pow(Bep,2)*pow(r,5)*pow(ssr,3)*pow(P,3) - 768*pow(Bep,2)*Bepp*pow(r,5)*pow(ssr,5)*pow(P,3) - 64*pow(Bep,2)*pow(r,2)*r_Der_Q*pow(ssr,2)*(-r + 8*Bep*Qr*r + 8*Bep*r*ssr*P) - 64*r_Der_P*(-3*pow(Bep,2)*pow(r,4)*pow(ssr,3) + 24*pow(Bep,3)*Qr*pow(r,4)*pow(ssr,3) + 2*pow(Bep,2)*pow(r,6)*pow(ssr,5) - 16*pow(Bep,3)*Qr*pow(r,6)*pow(ssr,5) + 20*pow(Bep,3)*pow(r,4)*pow(ssr,4)*P - 12*pow(Bep,3)*pow(r,6)*pow(ssr,6)*P) + 16*r_Der_ss*(Bep*pow(r,3)*ssr - 16*pow(Bep,2)*Qr*pow(r,3)*ssr + 64*pow(Bep,3)*pow(Qr,2)*pow(r,3)*ssr + 8*pow(Bep,2)*Qr*pow(r,5)*pow(ssr,3) - 64*pow(Bep,3)*pow(Qr,2)*pow(r,5)*pow(ssr,3) - 12*pow(Bep,2)*pow(r,3)*pow(ssr,2)*P + 96*pow(Bep,3)*Qr*pow(r,3)*pow(ssr,2)*P - 48*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4)*P + 32*pow(Bep,3)*pow(r,3)*pow(ssr,3)*pow(P,2))))/(pow(r,4) - 16*Bep*Qr*pow(r,4) + 64*pow(Bep,2)*pow(Qr,2)*pow(r,4) - 32*pow(Bep,2)*pow(r,4)*pow(ssr,4) + 256*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,6)*pow(ssr,4) + 256*pow(Bep,3)*pow(r,4)*r_Der_Q*pow(ssr,4) - 16*Bep*pow(r,4)*ssr*P + 128*pow(Bep,2)*Qr*pow(r,4)*ssr*P + 64*pow(Bep,2)*pow(r,4)*pow(ssr,2)*pow(P,2) + 128*r_Der_ss*(pow(Bep,2)*pow(r,4)*pow(ssr,3) - 8*pow(Bep,3)*Qr*pow(r,4)*pow(ssr,3) - 6*pow(Bep,3)*pow(r,4)*pow(ssr,4)*P) + (128*r_Der_nn*(pow(Bep,2)*pow(r,5)*pow(ssr,4) - 8*pow(Bep,3)*Qr*pow(r,5)*pow(ssr,4) - 6*pow(Bep,3)*pow(r,5)*pow(ssr,5)*P))/nn)
  ;
/*---------------------------------------------------------------------------*/
  double ep_dtp= (768*pow(Bep,3)*pow(r,4)*r_Der_Q*pow(ssr,4)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P) + (pow(r,4)*(-1 + 24*Bep*Qr - 192*pow(Bep,2)*pow(Qr,2) + 512*pow(Bep,3)*pow(Qr,3) - 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,2) + 256*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(ssr,2) + 96*pow(Bep,2)*pow(ssr,4) - 768*pow(Bep,3)*Qr*pow(ssr,4) - 768*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4) + 6144*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,4) + 28*Bep*ssr*P - 448*pow(Bep,2)*Qr*ssr*P + 1792*pow(Bep,3)*pow(Qr,2)*ssr*P + 192*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(ssr,3)*P - 768*pow(Bep,3)*pow(ssr,5)*P + 6144*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,5)*P - 288*pow(Bep,2)*pow(ssr,2)*pow(P,2) + 2304*pow(Bep,3)*Qr*pow(ssr,2)*pow(P,2) + 960*pow(Bep,3)*pow(ssr,3)*pow(P,3)))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P)
     ;
  double ep_dtq=
     0
  ;
  double ep_drp=
     -1536*pow(Bep,3)*pow(r,4)*r_Der_P*pow(ssr,4)*nn - (768*pow(Bep,3)*pow(r,5)*r_Der_Q*pow(ssr,5)*nn*(-1 + 4*Bep*Qr + 8*Bep*ssr*P))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P) - (pow(r,5)*ssr*nn*(1 - 36*Bep*Qr + 480*pow(Bep,2)*pow(Qr,2) - 2816*pow(Bep,3)*pow(Qr,3) + 6144*pow(Bep,4)*pow(Qr,4) + 32*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,2) - 256*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(ssr,2) - 96*pow(Bep,2)*pow(ssr,4) + 1152*pow(Bep,3)*Qr*pow(ssr,4) - 3072*pow(Bep,4)*pow(Qr,2)*pow(ssr,4) + 768*pow(Bep,2)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,4) - 9216*pow(Bep,3)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,4) + 24576*pow(Bep,4)*Bepp*pow(Qr,4)*pow(r,2)*pow(ssr,4) - 36*Bep*ssr*P + 1104*pow(Bep,2)*Qr*ssr*P - 10752*pow(Bep,3)*pow(Qr,2)*ssr*P + 33792*pow(Bep,4)*pow(Qr,3)*ssr*P + 1536*pow(Bep,2)*Bepp*Qr*pow(ssr,3)*P - 24576*pow(Bep,3)*Bepp*pow(Qr,2)*pow(ssr,3)*P + 98304*pow(Bep,4)*Bepp*pow(Qr,3)*pow(ssr,3)*P - 448*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(ssr,3)*P + 1280*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(ssr,3)*P + 1536*pow(Bep,3)*pow(ssr,5)*P - 9216*pow(Bep,4)*Qr*pow(ssr,5)*P - 12288*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,5)*P + 73728*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,5)*P + 512*pow(Bep,2)*pow(ssr,2)*pow(P,2) - 11520*pow(Bep,3)*Qr*pow(ssr,2)*pow(P,2) + 59392*pow(Bep,4)*pow(Qr,2)*pow(ssr,2)*pow(P,2) - 30720*pow(Bep,3)*Bepp*Qr*pow(ssr,4)*pow(P,2) + 245760*pow(Bep,4)*Bepp*pow(Qr,2)*pow(ssr,4)*pow(P,2) + 1536*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,4)*pow(P,2) - 6144*pow(Bep,4)*pow(ssr,6)*pow(P,2) + 49152*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,6)*pow(P,2) - 3264*pow(Bep,3)*pow(ssr,3)*pow(P,3) + 39168*pow(Bep,4)*Qr*pow(ssr,3)*pow(P,3) + 147456*pow(Bep,4)*Bepp*Qr*pow(ssr,5)*pow(P,3) + 7680*pow(Bep,4)*pow(ssr,4)*pow(P,4)))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*(-1 + 8*Bep*Qr + 12*Bep*ssr*P))
     ;
  double ep_drq=
     (768*pow(Bep,3)*pow(r,4)*pow(ssr,4)*t_Der_P*(-1 + 8*Bep*Qr + 8*Bep*ssr*P))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P) - (768*pow(r,5)*r_Der_P*(-(pow(Bep,3)*pow(ssr,5)*nn) + 4*pow(Bep,4)*Qr*pow(ssr,5)*nn + 8*pow(Bep,4)*pow(ssr,6)*nn*P))/(-1 + 8*Bep*Qr + 12*Bep*ssr*P) - (pow(r,4)*(nn - 32*Bep*Qr*nn + 384*pow(Bep,2)*pow(Qr,2)*nn - 2048*pow(Bep,3)*pow(Qr,3)*nn + 4096*pow(Bep,4)*pow(Qr,4)*nn - 48*pow(Bep,2)*pow(Qr,2)*pow(r,2)*pow(ssr,2)*nn + 768*pow(Bep,3)*pow(Qr,3)*pow(r,2)*pow(ssr,2)*nn - 3072*pow(Bep,4)*pow(Qr,4)*pow(r,2)*pow(ssr,2)*nn - 96*pow(Bep,2)*pow(ssr,4)*nn + 1536*pow(Bep,3)*Qr*pow(ssr,4)*nn - 6144*pow(Bep,4)*pow(Qr,2)*pow(ssr,4)*nn + 256*pow(Bep,4)*pow(Qr,4)*pow(r,4)*pow(ssr,4)*nn - 32*Bep*ssr*nn*P + 768*pow(Bep,2)*Qr*ssr*nn*P - 6144*pow(Bep,3)*pow(Qr,2)*ssr*nn*P + 16384*pow(Bep,4)*pow(Qr,3)*ssr*nn*P + 896*pow(Bep,3)*pow(Qr,2)*pow(r,2)*pow(ssr,3)*nn*P - 7168*pow(Bep,4)*pow(Qr,3)*pow(r,2)*pow(ssr,3)*nn*P + 1536*pow(Bep,3)*pow(ssr,5)*nn*P - 12288*pow(Bep,4)*Qr*pow(ssr,5)*nn*P + 3072*pow(Bep,3)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,5)*nn*P - 24576*pow(Bep,4)*Bepp*pow(Qr,3)*pow(r,2)*pow(ssr,5)*nn*P + 352*pow(Bep,2)*pow(ssr,2)*nn*pow(P,2) - 5632*pow(Bep,3)*Qr*pow(ssr,2)*nn*pow(P,2) + 22528*pow(Bep,4)*pow(Qr,2)*pow(ssr,2)*nn*pow(P,2) - 768*pow(Bep,2)*Bepp*pow(ssr,4)*nn*pow(P,2) + 12288*pow(Bep,3)*Bepp*Qr*pow(ssr,4)*nn*pow(P,2) - 49152*pow(Bep,4)*Bepp*pow(Qr,2)*pow(ssr,4)*nn*pow(P,2) - 3840*pow(Bep,4)*pow(Qr,2)*pow(r,2)*pow(ssr,4)*nn*pow(P,2) - 6144*pow(Bep,4)*pow(ssr,6)*nn*pow(P,2) - 24576*pow(Bep,4)*Bepp*pow(Qr,2)*pow(r,2)*pow(ssr,6)*nn*pow(P,2) - 1536*pow(Bep,3)*pow(ssr,3)*nn*pow(P,3) + 12288*pow(Bep,4)*Qr*pow(ssr,3)*nn*pow(P,3) + 12288*pow(Bep,3)*Bepp*pow(ssr,5)*nn*pow(P,3) - 98304*pow(Bep,4)*Bepp*Qr*pow(ssr,5)*nn*pow(P,3) + 2048*pow(Bep,4)*pow(ssr,4)*nn*pow(P,4) - 49152*pow(Bep,4)*Bepp*pow(ssr,6)*nn*pow(P,4)))/((-1 + 8*Bep*Qr + 8*Bep*ssr*P)*(-1 + 8*Bep*Qr + 12*Bep*ssr*P))
     ;
/*---------------------------------------------------------------------------*/
  double eq_dtp=
     0
  ;
  double eq_dtq=
     1
  ;
  double eq_drp=
  -   nn
  ;
  double eq_drq=
  -   nn*r*ssr
  ;
/*---------------------------------------------------------------------------*/
  double A= ep_dtp*eq_dtq- ep_dtq*eq_dtp;
  double B=
  -  (ep_dtp*eq_drq-ep_dtq*eq_drp)
  -  (ep_drp*eq_dtq-ep_drq*eq_dtp)
  ;
  double C= ep_drp*eq_drq-ep_drq*eq_drp;

  double D=  pow(B,2)-4*A*C;

  if (D<0) {
     ingoing_c=  0;
     outgoing_c= 0;
     return -1;
  }
/*---------------------------------------------------------------------------*/
  ingoing_c= (-B-sqrt(D))/(2*A);
  outgoing_c=(-B+sqrt(D))/(2*A);
/*---------------------------------------------------------------------------*/
  return 0;

}
//==============================================================================
void Diagnostics::check_for_elliptic_region(Grid_data &grid,
  const Field &n_v, const Field &s_v,
  const Field &p_v, const Field &q_v,
  const Field &phi_v, vector<double> &ingoing, vector<double> &outgoing)
  {

    vector<double> r = grid.r;
    int nx = grid.nx;
    vector<double> dr = grid.dr;
    double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;
    double ingoing_c=0., outgoing_c=0.;
    if(grid.exc_i==0){
   /*---------------------------------------------------------------------------*/
   /* r=0 */
      {
        int i = grid.exc_i;

        double nn = n_v.v[i];

        double r_Der_Q = Dx_ptpc_2nd(q_v.v[i+1], -q_v.v[i+1], dr[i]);
        double r_Der_ss = Dx_ptpc_2nd(s_v.v[i+1], -s_v.v[i+1], dr[i]);
        double P = p_v.v[i];

        double rr_Der_P =Dx_2_ptpc_2nd(p_v.v[i+1], p_v.v[i], p_v.v[i+1], dr[i]);
        double rr_Der_nn = Dx_2_ptpc_2nd(n_v.v[i+1], n_v.v[i], n_v.v[i+1], dr[i]);

        double Bep = beta_p(ls,lexp,mu, phi_v.v[i]);
        double Bepp = beta_pp(ls,lexp,mu, phi_v.v[i]);

        double ingoing_c=  0;
        double outgoing_c= 0;
        int status = 0.;

        double D = (-4*pow(nn,2)*(pow(-1 + 8*Bep*r_Der_Q,7) + 60*Bep*P*pow(1 - 8*Bep*r_Der_Q,6)*r_Der_ss - 1536*pow(Bep,2)*pow(P,2)*pow(1 - 8*Bep*r_Der_Q,5)*pow(r_Der_ss,2) + 21760*pow(Bep,3)*pow(P,3)*pow(1 - 8*Bep*r_Der_Q,4)*pow(r_Der_ss,3) + 16*pow(Bep,2)*pow(1 - 8*Bep*r_Der_Q,3)*(-11504*pow(Bep,2)*pow(P,4) + 48*pow(P,2)*(Bepp - 16*Bep*Bepp*r_Der_Q + 64*pow(Bep,2)*Bepp*pow(r_Der_Q,2)) + 4*(3 - 24*Bep*r_Der_Q - 192*pow(Bep,2)*pow(r_Der_Q,2) + 1536*pow(Bep,3)*pow(r_Der_Q,3)))*pow(r_Der_ss,4) - 128*pow(Bep,2)*pow(1 - 8*Bep*r_Der_Q,2)*(-7240*pow(Bep,3)*pow(P,5) + 264*Bep*pow(P,3)*(Bepp - 16*Bep*Bepp*r_Der_Q + 64*pow(Bep,2)*Bepp*pow(r_Der_Q,2)) + P*(54*Bep - 240*pow(Bep,2)*r_Der_Q - 6528*pow(Bep,3)*pow(r_Der_Q,2) + 39936*pow(Bep,4)*pow(r_Der_Q,3)))*pow(r_Der_ss,5) + 64*pow(Bep,2)*(1 - 8*Bep*r_Der_Q)*(3*(400*pow(Bep,2)*pow(P,2) + 3072*pow(Bep,2)*Bepp*pow(P,4) - 13312*pow(Bep,4)*pow(P,6)) - 96*Bep*(-84*pow(Bep,2)*pow(P,2) + 1536*pow(Bep,2)*Bepp*pow(P,4))*r_Der_Q + 64*pow(Bep,4)*P*(-5616*P + 9216*Bepp*pow(P,3))*pow(r_Der_Q,2) + 1744896*pow(Bep,5)*pow(P,2)*pow(r_Der_Q,3))*pow(r_Der_ss,6) - 256*pow(Bep,3)*(-11520*pow(Bep,4)*pow(P,7) + 19968*pow(Bep,2)*pow(P,5)*(Bepp - 16*Bep*Bepp*r_Der_Q + 64*pow(Bep,2)*Bepp*pow(r_Der_Q,2)) + 48*Bep*pow(P,3)*(-Bep + 1688*pow(Bep,2)*r_Der_Q - 26816*pow(Bep,3)*pow(r_Der_Q,2) + 107008*pow(Bep,4)*pow(r_Der_Q,3)))*pow(r_Der_ss,7) - 6144*pow(Bep,4)*(2*pow(Bep,2)*pow(P,6)*(-1792*Bepp + 14336*Bep*Bepp*r_Der_Q) - 4*Bep*pow(P,4)*(-252*Bep + 7616*pow(Bep,2)*r_Der_Q - 44800*pow(Bep,3)*pow(r_Der_Q,2)))*pow(r_Der_ss,8) - 786432*pow(Bep,4)*(-55*pow(Bep,3)*pow(P,5) + 48*pow(Bep,3)*Bepp*pow(P,7) + 632*pow(Bep,4)*pow(P,5)*r_Der_Q)*pow(r_Der_ss,9) - 94371840*pow(Bep,8)*pow(P,6)*pow(r_Der_ss,10)))/((1 - 8*Bep*r_Der_Q - 8*Bep*P*r_Der_ss)*pow(1 - 512*pow(Bep,3)*pow(r_Der_Q,3) + 288*pow(Bep,2)*pow(P,2)*pow(r_Der_ss,2) - 960*pow(Bep,3)*pow(P,3)*pow(r_Der_ss,3) + 96*pow(Bep,2)*(-1 + 8*Bep*r_Der_Q)*pow(r_Der_ss,4) + 32*pow(Bep,2)*pow(r_Der_Q,2)*(6 - 56*Bep*P*r_Der_ss) - 8*Bep*r_Der_Q*(3 - 56*Bep*P*r_Der_ss + 288*pow(Bep,2)*pow(P,2)*pow(r_Der_ss,2) + 96*pow(Bep,2)*(-1 + 8*Bep*r_Der_Q)*pow(r_Der_ss,4)) - 4*P*(7*Bep*r_Der_ss + 192*pow(Bep,3)*(-1 + 8*Bep*r_Der_Q)*pow(r_Der_ss,5)),2))
        ;

        double trace = 0.;

        if (D<0) {
          ingoing_c=  0;
          outgoing_c= 0;
          status = -1;
        }
        else{
          ingoing_c= 0.5*(trace - pow(D,0.5));
          outgoing_c=0.5*(trace + pow(D,0.5));
          status = 0.;
        }


        if (status==-1) {
          // cout<<"Elliptic region formation in flat space."<<endl;
          // cout<<"naked_elliptic_region at (i,r) = ("<<i<<","<<r[i]<<"), t = "<<grid.t_evolve<<endl;
          // std::exit(0);
        }
        ingoing[i]=   ingoing_c;
        outgoing[i]= outgoing_c;
      }
   /*---------------------------------------------------------------------------*/
   /* interior */
      for (int i=grid.exc_i + 1; i<nx-1; ++i) {
        double Bep=  beta_p(ls,lexp,mu,phi_v.v[i]);
        double Bepp= beta_pp(ls,lexp,mu,phi_v.v[i]);

        double r_Der_nn= Dx_ptpc_2nd(n_v.v[i+1], n_v.v[i-1], dr[i]);
        double r_Der_ss= Dx_ptpc_2nd(s_v.v[i+1], s_v.v[i-1], dr[i]);
        double r_Der_P= Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
        double r_Der_Q= Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);

        int status= compute_radial_characteristic(r[i],
        n_v.v[i], r_Der_nn,
        s_v.v[i], r_Der_ss,
        p_v.v[i], r_Der_P,
        q_v.v[i], r_Der_Q,
        Bep, Bepp,
        ingoing_c, outgoing_c);

        if (status==-1) {
          if(grid.NER==0){
          cout<<"Elliptic region formation in flat space."<<endl;
          cout<<"naked_elliptic_region at (i,r) = ("<<i<<","<<r[i]<<"), t = "<<grid.t_evolve<<endl;
          grid.NER = 1;
          grid.NER_index = i;
          }
          // std::exit(0);
        }
        ingoing[i]=   ingoing_c;
        outgoing[i]= outgoing_c;

      }
    }

    else{

  int new_exc_i = grid.exc_i;

  if(grid.bh_start==0){

  for(int i =grid.exc_i; i<grid.exc_i+1;i++) {


     double Bep=  beta_p(ls,lexp,mu,phi_v.v[i]);
     double Bepp= beta_pp(ls,lexp,mu,phi_v.v[i]);

     double r_Der_nn= Dx_ptp0_2nd(n_v.v[i+2], n_v.v[i+1], n_v.v[i], dr[i]);
     double r_Der_ss= Dx_ptp0_2nd(s_v.v[i+2], s_v.v[i+1], s_v.v[i], dr[i]);
     double r_Der_P= Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
     double r_Der_Q= Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);

     int status= compute_radial_characteristic(r[i],
     n_v.v[i], r_Der_nn,
     s_v.v[i], r_Der_ss,
     p_v.v[i], r_Der_P,
     q_v.v[i], r_Der_Q,
     Bep, Bepp,
     ingoing_c, outgoing_c);

     if (status==-1) {
       new_exc_i += 1;
     }
     ingoing[i]=   ingoing_c;
     outgoing[i]= outgoing_c;
   }

/*---------------------------------------------------------------------------*/
/* interior */
   for (int i=grid.exc_i+1; i<nx-1; ++i) {

     double Bep=  beta_p(ls,lexp,mu,phi_v.v[i]);
     double Bepp= beta_pp(ls,lexp,mu,phi_v.v[i]);

     double r_Der_nn= Dx_ptpc_2nd(n_v.v[i+1], n_v.v[i-1], dr[i]);
     double r_Der_ss= Dx_ptpc_2nd(s_v.v[i+1], s_v.v[i-1], dr[i]);
     double r_Der_P= Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
     double r_Der_Q= Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);


     int status= compute_radial_characteristic(r[i],
     n_v.v[i], r_Der_nn,
     s_v.v[i], r_Der_ss,
     p_v.v[i], r_Der_P,
     q_v.v[i], r_Der_Q,
     Bep, Bepp,
     ingoing_c, outgoing_c);

     if (status==-1) {
       new_exc_i += 1;
     }
     ingoing[i]=   ingoing_c;
     outgoing[i]= outgoing_c;

   }
 }
  else{

    for(int i =grid.exc_i; i<grid.ah_index-3;i++) {


       double Bep=  beta_p(ls,lexp,mu,phi_v.v[i]);
       double Bepp= beta_pp(ls,lexp,mu,phi_v.v[i]);

       double r_Der_nn= Dx_ptp0_2nd(n_v.v[i+2], n_v.v[i+1], n_v.v[i], dr[i]);
       double r_Der_ss= Dx_ptp0_2nd(s_v.v[i+2], s_v.v[i+1], s_v.v[i], dr[i]);
       double r_Der_P= Dx_ptp0_2nd(p_v.v[i+2], p_v.v[i+1], p_v.v[i], dr[i]);
       double r_Der_Q= Dx_ptp0_2nd(q_v.v[i+2], q_v.v[i+1], q_v.v[i], dr[i]);

       int status= compute_radial_characteristic(r[i],
       n_v.v[i], r_Der_nn,
       s_v.v[i], r_Der_ss,
       p_v.v[i], r_Der_P,
       q_v.v[i], r_Der_Q,
       Bep, Bepp,
       ingoing_c, outgoing_c);

       if (status==-1) {
         new_exc_i += 1;
       }
       ingoing[i]=   ingoing_c;
       outgoing[i]= outgoing_c;
     }

  /*---------------------------------------------------------------------------*/
  /* interior */
     for (int i=grid.ah_index-3; i<nx-1; ++i) {

       double Bep=  beta_p(ls,lexp,mu,phi_v.v[i]);
       double Bepp= beta_pp(ls,lexp,mu,phi_v.v[i]);

       double r_Der_nn= Dx_ptpc_2nd(n_v.v[i+1], n_v.v[i-1], dr[i]);
       double r_Der_ss= Dx_ptpc_2nd(s_v.v[i+1], s_v.v[i-1], dr[i]);
       double r_Der_P= Dx_ptpc_2nd(p_v.v[i+1], p_v.v[i-1], dr[i]);
       double r_Der_Q= Dx_ptpc_2nd(q_v.v[i+1], q_v.v[i-1], dr[i]);


       int status= compute_radial_characteristic(r[i],
       n_v.v[i], r_Der_nn,
       s_v.v[i], r_Der_ss,
       p_v.v[i], r_Der_P,
       q_v.v[i], r_Der_Q,
       Bep, Bepp,
       ingoing_c, outgoing_c);

       if (status==-1) {
         new_exc_i += 1;
       }
       ingoing[i]=   ingoing_c;
       outgoing[i]= outgoing_c;

     }

  }
  //=============================================================================
        if(new_exc_i>grid.exc_i)
        {
          grid.exc_i = new_exc_i;
          cout<<"Elliptic region at (i,r) = ("<<grid.exc_i<<","<<r[grid.exc_i]<<"), AH at (i,r) = ("<<grid.ah_index<<","<<r[grid.ah_index]<<"), t = "<<grid.t_evolve<<endl;
      }
        else{
          grid.exc_i = new_exc_i;
      }
      if (grid.exc_i >= grid.ah_index-3){
        // cout<<"Elliptic region within 3 grid points of ah; ";
        if(grid.NER==0){
        cout<<"Elliptic region at (i,r) = ("<<grid.exc_i<<","<<r[grid.exc_i]<<"), AH at (i,r) = ("<<grid.ah_index<<","<<r[grid.ah_index]<<"), t = "<<grid.t_evolve<<endl;
        cout<<"naked_elliptic_region outside horizon. "<<"Elliptic region at (i,r) = ("<<grid.exc_i<<","<<r[grid.exc_i]<<"), AH at (i,r) = ("<<grid.ah_index<<","<<r[grid.ah_index]<<"), t = "<<grid.t_evolve<<endl;
        grid.NER = 1;
        grid.NER_index = grid.exc_i;
        // std::exit(0);
      }
      }
  }

      return;


    }


//==============================================================================
double Diagnostics::e_rr_residual(double r, double nn, double r_Der_nn,
double ss, double r_Der_ss,
double P, double r_Der_P,
double Q, double r_Der_Q,
double Bep, double Bepp,
double t_Der_ss, double t_Der_P)
{

  double Qr = Q/r;
  double ssr = ss/r;
  double tm1 = (-2*r_Der_nn*(-1 + 8*Bep*Qr + 8*Bep*P*ssr))/(nn) - (2*t_Der_ss*(-1 + 8*Bep*Qr + 8*Bep*P*ssr))/(nn);
  double t0 = -pow(ssr,2) - (8*Bep*t_Der_P*pow(ssr,2))/nn - (pow(P,2)*(1 + 16*Bepp*pow(ssr,2)))/2. + r_Der_ss*(2*(-1 + 8*Bep*Qr)*ssr + 16*Bep*P*pow(ssr,2));
  double t1 = (8*Bep*r*r_Der_nn*Qr*pow(ssr,2))/nn + 8*Bep*r*r_Der_P*pow(ssr,3);
  double t2 = -0.5*(pow(r,2)*pow(Qr,2));

  return tm1/r + t0 + t1 + t2;
}
//==============================================================================
void Diagnostics::compute_e_rr_residual(Grid_data &grid, const vector<double> &n_v, const vector<double> &s_v,
  const vector<double> &p_v, const vector<double> &q_v, const vector<double> &phi_v,
  const vector<double> &n_v_np1,
  const vector<double> &s_v_np1,
  const vector<double> &p_v_np1,
  const vector<double> &q_v_np1,const vector<double> &phi_v_np1, vector<double> &residual)
  {

    vector<double> dr = grid.dr;
    vector<double> r = grid.r;
    double dt = grid.dt;
    double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;
    int nx = grid.nx;


    if(grid.exc_i==0){

      {
        int i = grid.exc_i;

        double rr_Der_nn= Dx_2_ptpc_2nd((n_v[i+1] + n_v_np1[i+1])/2., (n_v[i] + n_v_np1[i])/2., (n_v[i+1] + n_v_np1[i+1])/2., dr[i]);

        double nn = (n_v[i] + n_v_np1[i])/2.;

        double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], -s_v[i+1], dr[i])
                    + Dx_ptpc_2nd(s_v_np1[i+1], -s_v_np1[i+1], dr[i]))/2.;
        double P = (p_v[i] + p_v_np1[i])/2.;
        double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], -q_v[i+1], dr[i])
                    + Dx_ptpc_2nd(q_v_np1[i+1], -q_v_np1[i+1], dr[i]))/2.;

        double s1dot = (Dx_ptpc_2nd(s_v_np1[i+1], -s_v_np1[i+1], dr[i]) - Dx_ptpc_2nd(s_v[i+1], -s_v[i+1], dr[i]))/dt;


        double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

        residual[i] = -0.5*(4*rr_Der_nn*(-1 + 8*Bep*r_Der_Q + 8*Bep*P*r_Der_ss) + nn*(2*(3 - 16*Bep*r_Der_Q)*pow(r_Der_ss,2) - 32*Bep*P*pow(r_Der_ss,3) + pow(P,2)*(1 + 16*Bepp*pow(r_Der_ss,2))) + 4*((-1 + 8*Bep*r_Der_Q)*s1dot + 8*Bep*P*r_Der_ss*s1dot + 4*Bep*pow(r_Der_ss,2)*t_Der_P))/nn;

      }

      for(int i = grid.exc_i+1; i<nx-1; i++){

        double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                    + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

        double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
                    + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;

        double r_Der_P= (Dx_ptpc_2nd(p_v[i+1], p_v[i-1], dr[i])
                    + Dx_ptpc_2nd(p_v_np1[i+1], p_v_np1[i-1], dr[i]))/2.;

        double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], q_v[i-1], dr[i])
                    + Dx_ptpc_2nd(q_v_np1[i+1], q_v_np1[i-1], dr[i]))/2.;

        double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

        residual[i] = e_rr_residual(r[i], (n_v[i] + n_v_np1[i])/2.,
        r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
        (p_v[i] + p_v_np1[i])/2.,
        r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
      }
      residual[nx-1] = 0.;

    }

    else{
      if(grid.bh_start==0){
      {
        int i = grid.exc_i;

        double r_Der_nn= (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i])
                    + Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2. ;

        double r_Der_ss= (Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i])
                    + Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i]))/2.;

        double r_Der_P= (Dx_ptp0_2nd(p_v[i+2], p_v[i+1], p_v[i], dr[i])
                    + Dx_ptp0_2nd(p_v_np1[i+2], p_v_np1[i+1], p_v_np1[i], dr[i]))/2.;

        double r_Der_Q= (Dx_ptp0_2nd(q_v[i+2], q_v[i+1], q_v[i], dr[i])
                    + Dx_ptp0_2nd(q_v_np1[i+2], q_v_np1[i+1], q_v_np1[i], dr[i]))/2.;

        double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

        residual[i] = e_rr_residual(r[i], (n_v[i] + n_v_np1[i])/2.,
        r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
        (p_v[i] + p_v_np1[i])/2.,
        r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
      }

      for(int i = grid.exc_i+1; i<nx-1; i++){

        double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                    + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

        double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
                    + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;

        double r_Der_P= (Dx_ptpc_2nd(p_v[i+1], p_v[i-1], dr[i])
                    + Dx_ptpc_2nd(p_v_np1[i+1], p_v_np1[i-1], dr[i]))/2.;

        double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], q_v[i-1], dr[i])
                    + Dx_ptpc_2nd(q_v_np1[i+1], q_v_np1[i-1], dr[i]))/2.;

        double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

        residual[i] = e_rr_residual(r[i], (n_v[i] + n_v_np1[i])/2.,
        r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
        (p_v[i] + p_v_np1[i])/2.,
        r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
      }
      residual[nx-1] = 0.;
    }
      else{
      for(int i = grid.exc_i; i<grid.ah_index-3; i++){

        double r_Der_nn= (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i])
                    + Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2. ;

        double r_Der_ss= (Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i])
                    + Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i]))/2.;

        double r_Der_P= (Dx_ptp0_2nd(p_v[i+2], p_v[i+1], p_v[i], dr[i])
                    + Dx_ptp0_2nd(p_v_np1[i+2], p_v_np1[i+1], p_v_np1[i], dr[i]))/2.;

        double r_Der_Q= (Dx_ptp0_2nd(q_v[i+2], q_v[i+1], q_v[i], dr[i])
                    + Dx_ptp0_2nd(q_v_np1[i+2], q_v_np1[i+1], q_v_np1[i], dr[i]))/2.;

        double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

        residual[i] = e_rr_residual(r[i], (n_v[i] + n_v_np1[i])/2.,
        r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
        (p_v[i] + p_v_np1[i])/2.,
        r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
        }

      for(int i = grid.ah_index-3; i<nx-1; i++){

        double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                    + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

        double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
                    + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;

        double r_Der_P= (Dx_ptpc_2nd(p_v[i+1], p_v[i-1], dr[i])
                    + Dx_ptpc_2nd(p_v_np1[i+1], p_v_np1[i-1], dr[i]))/2.;

        double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], q_v[i-1], dr[i])
                    + Dx_ptpc_2nd(q_v_np1[i+1], q_v_np1[i-1], dr[i]))/2.;

        double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

        double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

        residual[i] = e_rr_residual(r[i], (n_v[i] + n_v_np1[i])/2.,
        r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
        (p_v[i] + p_v_np1[i])/2.,
        r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      }
      residual[nx-1] = 0.;
    }


    }



  }

//==============================================================================
double Diagnostics::get_NCC_in(double nn, double r_Der_nn, double ss, double t_Der_ss)
{
  return t_Der_ss*nn + r_Der_nn*nn*pow(1 + ss,2);
}
//==============================================================================
double Diagnostics::get_NCC_out(double nn, double r_Der_nn, double ss, double t_Der_ss)
{
  return t_Der_ss*nn + r_Der_nn*nn*pow(-1 + ss,2);
}
//==============================================================================
void Diagnostics::compute_NCC(Grid_data &grid,
  const std::vector<double> &n_v,
  const std::vector<double> &s_v,
  const std::vector<double> &n_v_np1,
  const std::vector<double> &s_v_np1,
  std::vector<double> &ncc_in,
  std::vector<double> &ncc_out)
  {
    vector<double> dr = grid.dr;
    // vector<double> r = grid.r;
    double dt = grid.dt;
    int nx = grid.nx;
    double r_Der_nn = 0., t_Der_ss = 0.;

    if(grid.exc_i==0){

      {
        int i = grid.exc_i;

        r_Der_nn=0. ;

        t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        ncc_in[i] = get_NCC_in((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);
        ncc_out[i] = get_NCC_out((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);

      }

      for(int i = grid.exc_i+1; i<nx-1; i++){

        r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                    + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

        t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        ncc_in[i] = get_NCC_in((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);
        ncc_out[i] = get_NCC_out((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);

      }
      ncc_in[nx-1] = 0.;
      ncc_out[nx-1] = 0.;

    }

    else{

      if(grid.bh_start==0){

      {
        int i = grid.exc_i;

        r_Der_nn= (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i]) +
                    Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2.;

        t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        ncc_in[i] = get_NCC_in((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);
        ncc_out[i] = get_NCC_out((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);

      }

      for(int i = grid.exc_i+1; i<nx-1; i++){

        r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                    + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

        t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

        ncc_in[i] = get_NCC_in((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);
        ncc_out[i] = get_NCC_out((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);

      }
      ncc_in[nx-1] = 0.;
      ncc_out[nx-1] = 0.;



    }
      else{


          for(int i = grid.exc_i; i< grid.ah_index-3; i++){


            r_Der_nn= (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i]) +
                        Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2.;

            t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

            ncc_in[i] = get_NCC_in((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);
            ncc_out[i] = get_NCC_out((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);

          }

          for(int i = grid.ah_index-3; i<nx-1; i++){

            r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                        + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

            t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

            ncc_in[i] = get_NCC_in((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);
            ncc_out[i] = get_NCC_out((n_v[i] +n_v_np1[i])/2., r_Der_nn, (s_v[i] +s_v_np1[i])/2., t_Der_ss);

          }
          ncc_in[nx-1] = 0.;
          ncc_out[nx-1] = 0.;

  }


  }
}
//==============================================================================
double Diagnostics::get_GB_Val_reduced(double r, double nn, double r_Der_nn, double ss, double r_Der_ss, double P, double r_Der_P, double Q, double r_Der_Q, double Bep, double Bepp, double t_Der_ss, double t_Der_P)
{
  double ssr = ss/r;
  double Qr = Q/r;
  return
  (64*pow(r_Der_nn,2)*pow(ssr,3)*(Bep*Qr*pow(r,2)*ssr + Bep*P))/(pow(nn,2)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (16*pow(r_Der_ss,2)*pow(ssr,2)*(-1 + 8*Bep*Qr + 4*Bep*ssr*P))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) - (4*pow(ssr,2)*(pow(Qr,2)*pow(r,2) - pow(P,2)))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) + t_Der_P*((64*Bep*r*r_Der_nn*pow(ssr,4))/(pow(nn,2)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (64*Bep*r_Der_ss*pow(ssr,3))/(nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P))) + t_Der_ss*((8*(-1 + 8*Bepp*pow(Qr,2)*pow(r,2))*pow(ssr,2))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (64*Bep*r_Der_Q*pow(ssr,2))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) - (16*r_Der_nn*pow(ssr,2)*(-1 + 8*Bep*Qr + 4*Bep*ssr*P))/(pow(nn,2)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) - (16*r_Der_ss*ssr*(-1 + 8*Bep*Qr + 4*Bep*ssr*P))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P))) + r_Der_ss*((-64*Bep*r_Der_Q*pow(ssr,3))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) - (64*Bep*r*r_Der_P*pow(ssr,4))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) - (16*pow(ssr,3)*(-1 + 4*Bepp*pow(Qr,2)*pow(r,2) - 4*Bepp*pow(P,2)))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + r_Der_nn*((64*Bep*r_Der_Q*pow(ssr,2))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) - (64*r_Der_P*pow(ssr,3)*(-2*Bep + Bep*pow(r,2)*pow(ssr,2)))/(nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (16*r_Der_ss*ssr*(1 - 8*Bep*Qr - pow(r,2)*pow(ssr,2) + 12*Bep*Qr*pow(r,2)*pow(ssr,2) - 4*Bep*ssr*P + 4*Bep*pow(r,2)*pow(ssr,3)*P))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (8*pow(ssr,2)*(-1 + 8*Bepp*pow(Qr,2)*pow(r,2) + pow(r,2)*pow(ssr,2) + 16*Bepp*Qr*pow(r,2)*ssr*P + 8*Bepp*pow(r,2)*pow(ssr,2)*pow(P,2)))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)))
  ;
}
//==============================================================================
double Diagnostics::get_Ricci_Val_reduced(double r, double nn, double r_Der_nn, double ss, double r_Der_ss, double P, double r_Der_P, double Q, double r_Der_Q, double Bep, double Bepp, double t_Der_ss, double t_Der_P)
{
  double ssr = ss/r;
  double Qr = Q/r;
  return
  (-16*Bep*pow(r_Der_ss,2)*ssr*P)/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) + (16*pow(r_Der_nn,2)*(Bep*Qr*pow(r,2)*pow(ssr,2) + Bep*ssr*P))/(pow(nn,2)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (-(pow(Qr,2)*pow(r,2)) - 2*pow(ssr,2) + 16*Bep*Qr*pow(ssr,2) + 16*Bep*pow(ssr,3)*P + pow(P,2))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) + t_Der_P*((16*Bep*r*r_Der_nn*pow(ssr,2))/(pow(nn,2)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (16*Bep*r_Der_ss*ssr)/(nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P))) + t_Der_ss*((16*Bep*r_Der_Q)/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (16*Bep*r_Der_nn*ssr*P)/(pow(nn,2)*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (16*Bep*r_Der_ss*P)/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (2*(1 - 16*Bep*Qr + 8*Bepp*pow(Qr,2)*pow(r,2) - 16*Bep*ssr*P))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P))) + r_Der_ss*((-16*Bep*r_Der_Q*ssr)/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) - (16*Bep*r*r_Der_P*pow(ssr,2))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P) - (4*(ssr - 16*Bep*Qr*ssr + 4*Bepp*pow(Qr,2)*pow(r,2)*ssr - 16*Bep*pow(ssr,2)*P - 4*Bepp*ssr*pow(P,2)))/(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + r_Der_nn*((16*Bep*r_Der_Q)/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) - (16*r_Der_P*(-2*Bep*ssr + Bep*pow(r,2)*pow(ssr,3)))/(nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (16*r_Der_ss*(Bep*Qr*pow(r,2)*ssr + Bep*P - Bep*pow(r,2)*pow(ssr,2)*P))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)) + (2*(1 - 16*Bep*Qr + 8*Bepp*pow(Qr,2)*pow(r,2) - pow(r,2)*pow(ssr,2) + 16*Bep*Qr*pow(r,2)*pow(ssr,2) - 16*Bep*ssr*P + 16*Bepp*Qr*pow(r,2)*ssr*P + 16*Bep*pow(r,2)*pow(ssr,3)*P + 8*Bepp*pow(r,2)*pow(ssr,2)*pow(P,2)))/(r*nn*(-1 + 8*Bep*Qr + 8*Bep*ssr*P)))
  ;
}
//==============================================================================
void Diagnostics::compute_GB_Ricci(Grid_data &grid, const vector<double> &n_v, const vector<double> &s_v,
  const vector<double> &p_v, const vector<double> &q_v, const vector<double> &phi_v,
  const vector<double> &n_v_np1,
  const vector<double> &s_v_np1,
  const vector<double> &p_v_np1,
  const vector<double> &q_v_np1,const vector<double> &phi_v_np1,
  vector<double> &gb,
  vector<double> &ricci)
{

  vector<double> dr = grid.dr;
  vector<double> r = grid.r;
  double dt = grid.dt;
  int nx = grid.nx;
  int exc_i = grid.exc_i;
  double ls = grid.ls, lexp = grid.lexp, mu = grid.mu;

  if(exc_i ==0){

    {
      int i = exc_i;

      double nn = (n_v[i] + n_v_np1[i])/2.;

      double r_Der_nn = 0.;

      double rr_Der_nn = (Dx_2_ptpc_2nd(n_v_np1[i+1], n_v_np1[i], n_v_np1[i+1], dr[i]) +
                          Dx_2_ptpc_2nd(n_v[i+1], n_v[i], n_v[i+1], dr[i]))/2.;

      double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;

      double tr_Der_nn = 0.;

      double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], -s_v[i+1], dr[i])
                  + Dx_ptpc_2nd(s_v_np1[i+1], -s_v_np1[i+1], dr[i]))/2.;

      double rr_Der_ss = (Dx_2_ptpc_2nd(s_v_np1[i+1], s_v_np1[i], -s_v_np1[i+1], dr[i]) +
                          Dx_2_ptpc_2nd(s_v[i+1], s_v[i], -s_v[i+1], dr[i]))/2.;

      double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

      double tr_Der_ss = (Dx_ptpc_2nd(s_v_np1[i+1], -s_v_np1[i+1], dr[i])
                          - Dx_ptpc_2nd(s_v[i+1], -s_v[i+1], dr[i]))/dt;

      gb[i] = (24*(-(rr_Der_nn*pow(r_Der_ss,2)) + nn*pow(r_Der_ss,4) - pow(r_Der_ss,2)*tr_Der_ss))/nn;

      ricci[i] = (6*(-rr_Der_nn + 2*nn*pow(r_Der_ss,2) - tr_Der_ss))/nn;
    }

    for(int i = exc_i + 1; i<nx-1; i++){
      double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                  + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

      double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
                  + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;

      double r_Der_P= (Dx_ptpc_2nd(p_v[i+1], p_v[i-1], dr[i])
                  + Dx_ptpc_2nd(p_v_np1[i+1], p_v_np1[i-1], dr[i]))/2.;

      double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], q_v[i-1], dr[i])
                  + Dx_ptpc_2nd(q_v_np1[i+1], q_v_np1[i-1], dr[i]))/2.;

      double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

      double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

      gb[i] = get_GB_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      ricci[i] = get_Ricci_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      //cout<<r[i]<<"\t"<<ricci[i]<<"\t"<<r_Der_nn<<"\t"<<rr_Der_nn<<"\t"<<t_Der_nn<<"\t"<<tr_Der_nn<<"\t"<<r_Der_ss<<"\t"<<rr_Der_ss<<"\t"<<t_Der_ss<<"\t"<<tr_Der_ss<<endl;
    }
    { int i = nx-1;

      gb[i] = 0.;
      ricci[i] = 0.;

    }
  }


  else{
    if(grid.bh_start==0){
    {
      int i = grid.exc_i;

      double r_Der_nn= (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i])
                  + Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2. ;

      double r_Der_ss= (Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i])
                  + Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i]))/2.;

      double r_Der_P= (Dx_ptp0_2nd(p_v[i+2], p_v[i+1], p_v[i], dr[i])
                  + Dx_ptp0_2nd(p_v_np1[i+2], p_v_np1[i+1], p_v_np1[i], dr[i]))/2.;

      double r_Der_Q= (Dx_ptp0_2nd(q_v[i+2], q_v[i+1], q_v[i], dr[i])
                  + Dx_ptp0_2nd(q_v_np1[i+2], q_v_np1[i+1], q_v_np1[i], dr[i]))/2.;

      double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

      double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

      gb[i] = get_GB_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      ricci[i] = get_Ricci_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
    }

    for(int i = grid.exc_i+1; i<nx-1; i++){

      double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                  + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

      double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
                  + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;

      double r_Der_P= (Dx_ptpc_2nd(p_v[i+1], p_v[i-1], dr[i])
                  + Dx_ptpc_2nd(p_v_np1[i+1], p_v_np1[i-1], dr[i]))/2.;

      double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], q_v[i-1], dr[i])
                  + Dx_ptpc_2nd(q_v_np1[i+1], q_v_np1[i-1], dr[i]))/2.;

      double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

      double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

      gb[i] = get_GB_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      ricci[i] = get_Ricci_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
    }
    { int i = nx-1;

      gb[i] = 0.;
      ricci[i] = 0.;

    }
  }
    else{
    for(int i = grid.exc_i; i<grid.ah_index-3; i++){

      double r_Der_nn= (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i])
                  + Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2. ;

      double r_Der_ss= (Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i])
                  + Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i]))/2.;

      double r_Der_P= (Dx_ptp0_2nd(p_v[i+2], p_v[i+1], p_v[i], dr[i])
                  + Dx_ptp0_2nd(p_v_np1[i+2], p_v_np1[i+1], p_v_np1[i], dr[i]))/2.;

      double r_Der_Q= (Dx_ptp0_2nd(q_v[i+2], q_v[i+1], q_v[i], dr[i])
                  + Dx_ptp0_2nd(q_v_np1[i+2], q_v_np1[i+1], q_v_np1[i], dr[i]))/2.;

      double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

      double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

      gb[i] = get_GB_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      ricci[i] = get_Ricci_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);
      }

    for(int i = grid.ah_index-3; i<nx-1; i++){

      double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
                  + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;

      double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
                  + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;

      double r_Der_P= (Dx_ptpc_2nd(p_v[i+1], p_v[i-1], dr[i])
                  + Dx_ptpc_2nd(p_v_np1[i+1], p_v_np1[i-1], dr[i]))/2.;

      double r_Der_Q= (Dx_ptpc_2nd(q_v[i+1], q_v[i-1], dr[i])
                  + Dx_ptpc_2nd(q_v_np1[i+1], q_v_np1[i-1], dr[i]))/2.;

      double Bep = beta_p(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double Bepp = beta_pp(ls,lexp,mu, (phi_v[i] +phi_v_np1[i])/2. );

      double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;

      double t_Der_P = (p_v_np1[i] - p_v[i])/dt;

      gb[i] = get_GB_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

      ricci[i] = get_Ricci_Val_reduced(r[i], (n_v[i] + n_v_np1[i])/2.,
      r_Der_nn, (s_v[i] + s_v_np1[i])/2., r_Der_ss,
      (p_v[i] + p_v_np1[i])/2.,
      r_Der_P, (q_v[i] + q_v_np1[i])/2., r_Der_Q, Bep, Bepp, t_Der_ss, t_Der_P);

    }
    { int i = nx-1;

      gb[i] = 0.;
      ricci[i] = 0.;

    }
  }
}
}
//============================================================================================
double Diagnostics::get_GB_Val(double r, double nn, double r_Der_nn, double rr_Der_nn, double t_Der_nn, double tr_Der_nn, double ss, double r_Der_ss, double rr_Der_ss, double t_Der_ss, double tr_Der_ss)
{
  double ans = 0;
  double ssr = ss/r;
  // ans = 8*pow(r,3)*r_Der_nn*pow(ssr,3)*t_Der_nn - 8*pow(r,3)*pow(ssr,3)*tr_Der_nn*nn + 8*r*r_Der_nn*r_Der_ss*ssr*(-2 + 5*pow(r,2)*pow(ssr,2))*pow(nn,2) - 8*pow(r,2)*pow(ssr,2)*tr_Der_ss*pow(nn,2) + 24*pow(r,2)*pow(r_Der_ss,2)*pow(ssr,2)*pow(nn,3) + 8*pow(r,3)*rr_Der_ss*pow(ssr,3)*pow(nn,3) + t_Der_ss*(-24*pow(r,2)*r_Der_nn*pow(ssr,2)*nn - 16*r*r_Der_ss*ssr*pow(nn,2)) + 8*r*rr_Der_nn*ssr*(-(r*ssr*pow(nn,2)) + pow(r,3)*pow(ssr,3)*pow(nn,2))
  // ;
  ans = 8*pow(r,3)*r_Der_nn*pow(ssr,3)*t_Der_nn - 8*pow(r,3)*pow(ssr,3)*tr_Der_nn*nn - 8*pow(r,2)*pow(ssr,2)*tr_Der_ss*pow(nn,2) + 24*pow(r,2)*pow(r_Der_ss,2)*pow(ssr,2)*pow(nn,3) + 8*pow(r,3)*rr_Der_ss*pow(ssr,3)*pow(nn,3) + t_Der_ss*(-24*pow(r,2)*r_Der_nn*pow(ssr,2)*nn - 16*r*r_Der_ss*ssr*pow(nn,2)) + 8*r*r_Der_nn*r_Der_ss*ssr*(-2*pow(nn,2) + 5*pow(r,2)*pow(ssr,2)*pow(nn,2)) + 8*r*rr_Der_nn*ssr*(-(r*ssr*pow(nn,2)) + pow(r,3)*pow(ssr,3)*pow(nn,2))
  ;
  ans /= pow(r,2)*pow(nn,3);

  return ans;
}
// //==============================================================================
double Diagnostics::get_Ricci_Val(double r, double nn, double r_Der_nn, double rr_Der_nn, double t_Der_nn, double tr_Der_nn, double ss, double r_Der_ss, double rr_Der_ss, double t_Der_ss, double tr_Der_ss)
{
  double ans = 0.;
  double ssr = ss/r;

  ans = 2*pow(r,3)*r_Der_nn*ssr*t_Der_nn - 2*pow(r,3)*ssr*tr_Der_nn*nn - 2*pow(r,2)*tr_Der_ss*pow(nn,2) + 2*pow(r,2)*pow(r_Der_ss,2)*pow(nn,3) + 2*pow(r,3)*rr_Der_ss*ssr*pow(nn,3) + 8*pow(r,2)*r_Der_ss*ssr*pow(nn,3) + 2*pow(r,2)*pow(ssr,2)*pow(nn,3) + t_Der_ss*(-2*pow(r,2)*r_Der_nn*nn - 4*r*pow(nn,2)) + 2*rr_Der_nn*(-(pow(r,2)*pow(nn,2)) + pow(r,4)*pow(ssr,2)*pow(nn,2)) + r_Der_nn*(6*pow(r,3)*r_Der_ss*ssr*pow(nn,2) + 4*(-(r*pow(nn,2)) + pow(r,3)*pow(ssr,2)*pow(nn,2)))
  ;

  ans /= pow(r,2)*pow(nn,3);

  return ans;
}
// //==============================================================================
// void Diagnostics::compute_GB_Ricci(Grid_data &grid,
// const vector<double> &n_v,
// const vector<double> &s_v,
// const vector<double> &n_v_np1,
// const vector<double> &s_v_np1,
// vector<double> &gb,
// vector<double> &ricci)
// {
//
//   vector<double> dr = grid.dr;
//   vector<double> r = grid.r;
//   double dt = grid.dt;
//   int nx = grid.nx;
//   int exc_i = grid.exc_i;
//
//   if(exc_i ==0){
//
//     {
//       int i = exc_i;
//
//       double nn = (n_v[i] + n_v_np1[i])/2.;
//
//       double r_Der_nn = 0.;
//
//       double rr_Der_nn = (Dx_2_ptpc_2nd(n_v_np1[i+1], n_v_np1[i], n_v_np1[i+1], dr[i]) +
//                           Dx_2_ptpc_2nd(n_v[i+1], n_v[i], n_v[i+1], dr[i]))/2.;
//
//       double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;
//
//       double tr_Der_nn = 0.;
//
//       double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], -s_v[i+1], dr[i])
//                   + Dx_ptpc_2nd(s_v_np1[i+1], -s_v_np1[i+1], dr[i]))/2.;
//
//       double rr_Der_ss = (Dx_2_ptpc_2nd(s_v_np1[i+1], s_v_np1[i], -s_v_np1[i+1], dr[i]) +
//                           Dx_2_ptpc_2nd(s_v[i+1], s_v[i], -s_v[i+1], dr[i]))/2.;
//
//       double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;
//
//       double tr_Der_ss = (Dx_ptpc_2nd(s_v_np1[i+1], -s_v_np1[i+1], dr[i])
//                           - Dx_ptpc_2nd(s_v[i+1], -s_v[i+1], dr[i]))/dt;
//
//       gb[i] = (24*(-(rr_Der_nn*pow(r_Der_ss,2)) + nn*pow(r_Der_ss,4) - pow(r_Der_ss,2)*tr_Der_ss))/nn;
//
//       ricci[i] = (6*(-rr_Der_nn + 2*nn*pow(r_Der_ss,2) - tr_Der_ss))/nn;
//     }
//
//     for(int i = exc_i + 1; i<nx-1; i++){
//       double r_Der_nn= (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
//                   + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2. ;
//
//       double rr_Der_nn = (Dx_2_ptpc_2nd(n_v_np1[i+1], n_v_np1[i], n_v_np1[i-1], dr[i]) +
//                           Dx_2_ptpc_2nd(n_v[i+1], n_v[i], n_v[i-1], dr[i]))/2.;
//
//       double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;
//
//       double tr_Der_nn = (Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i])
//                           - Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i]))/dt;
//
//       double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
//                   + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;
//
//       double rr_Der_ss = (Dx_2_ptpc_2nd(s_v_np1[i+1], s_v_np1[i], s_v_np1[i-1], dr[i]) +
//                           Dx_2_ptpc_2nd(s_v[i+1], s_v[i], s_v[i-1], dr[i]))/2.;
//
//       double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;
//
//       double tr_Der_ss = (Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i])
//                           - Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i]))/dt;
//
//       gb[i] = get_GB_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//         t_Der_nn, tr_Der_nn,
//         (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//
//       ricci[i] = get_Ricci_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//         t_Der_nn, tr_Der_nn,
//         (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//
//       //cout<<r[i]<<"\t"<<ricci[i]<<"\t"<<r_Der_nn<<"\t"<<rr_Der_nn<<"\t"<<t_Der_nn<<"\t"<<tr_Der_nn<<"\t"<<r_Der_ss<<"\t"<<rr_Der_ss<<"\t"<<t_Der_ss<<"\t"<<tr_Der_ss<<endl;
//     }
//     { int i = nx-1;
//
//       gb[i] = 0.;
//       ricci[i] = 0.;
//
//     }
//   }
//
//
//   else{
//     if(grid.bh_start ==0){
//       {
//         int i = exc_i;
//         double r_Der_nn = (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i])
//                     + Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2. ;
//
//         double rr_Der_nn = 0.5*(Dx_2_ptp0_2nd(n_v[i+3],n_v[i+2],n_v[i+1], n_v[i], dr[i]) + Dx_2_ptp0_2nd(n_v_np1[i+3],n_v_np1[i+2],n_v_np1[i+1], n_v_np1[i], dr[i]) );
//
//         double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;
//
//         double tr_Der_nn = ( Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]) - Dx_2_ptp0_2nd(n_v_np1[i+3],n_v_np1[i+2],n_v_np1[i+1], n_v_np1[i], dr[i]) )/dt ;
//
//         double r_Der_ss= (Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i])
//                     + Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i]))/2.;
//
//         double rr_Der_ss = 0.5*(Dx_2_ptp0_2nd(s_v[i+3],s_v[i+2],s_v[i+1], s_v[i], dr[i]) + Dx_2_ptp0_2nd(s_v_np1[i+3],s_v_np1[i+2],s_v_np1[i+1], s_v_np1[i], dr[i]) );
//
//         double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;
//
//         double tr_Der_ss = (Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i])
//                             - Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i]) )/dt;
//
//         gb[i] = get_GB_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//
//         ricci[i] = get_Ricci_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//       }
//       for(int i = exc_i + 1; i<nx-1; i++){
//         double r_Der_nn = (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
//                     + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2.;
//
//         double rr_Der_nn = 0.5*(Dx_2_ptpc_2nd(n_v[i+1], n_v[i], n_v[i-1], dr[i]) + Dx_2_ptpc_2nd(n_v_np1[i+1], n_v_np1[i], n_v_np1[i-1], dr[i]) );
//
//         double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;
//
//         double tr_Der_nn = (Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i])
//                             - Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i]))/dt;;
//
//         double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
//                     + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;
//
//         double rr_Der_ss = 0.5*(Dx_2_ptpc_2nd(s_v[i+1], s_v[i], s_v[i-1], dr[i]) + Dx_2_ptpc_2nd(s_v_np1[i+1], s_v_np1[i], s_v_np1[i-1], dr[i]) );
//
//         double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;
//
//         double tr_Der_ss = (Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i])
//                             - Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i]))/dt;
//
//         gb[i] = get_GB_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//
//         ricci[i] = get_Ricci_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//       }
//       { int i = nx-1;
//
//         gb[i] = 0.;
//         ricci[i] = 0.;
//
//       }
//
//
//     }
//     else{
//
//       for(int i = grid.exc_i; i<grid.ah_index-3; i++){
//
//         double r_Der_nn = (Dx_ptp0_2nd(n_v[i+2], n_v[i+1], n_v[i], dr[i])
//                     + Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]))/2. ;
//
//         double rr_Der_nn = 0.5*(Dx_2_ptp0_2nd(n_v[i+3],n_v[i+2],n_v[i+1], n_v[i], dr[i]) + Dx_2_ptp0_2nd(n_v_np1[i+3],n_v_np1[i+2],n_v_np1[i+1], n_v_np1[i], dr[i]) );
//
//         double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;
//
//         double tr_Der_nn = ( Dx_ptp0_2nd(n_v_np1[i+2], n_v_np1[i+1], n_v_np1[i], dr[i]) - Dx_2_ptp0_2nd(n_v_np1[i+3],n_v_np1[i+2],n_v_np1[i+1], n_v_np1[i], dr[i]) )/dt ;
//
//         double r_Der_ss= (Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i])
//                     + Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i]))/2.;
//
//         double rr_Der_ss = 0.5*(Dx_2_ptp0_2nd(s_v[i+3],s_v[i+2],s_v[i+1], s_v[i], dr[i]) + Dx_2_ptp0_2nd(s_v_np1[i+3],s_v_np1[i+2],s_v_np1[i+1], s_v_np1[i], dr[i]) );
//
//         double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;
//
//         double tr_Der_ss = (Dx_ptp0_2nd(s_v_np1[i+2], s_v_np1[i+1], s_v_np1[i], dr[i])
//                             - Dx_ptp0_2nd(s_v[i+2], s_v[i+1], s_v[i], dr[i]) )/dt;
//
//         gb[i] = get_GB_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//
//         ricci[i] = get_Ricci_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//       }
//
//       for(int i = grid.ah_index-3; i<nx-1; i++){
//
//         double r_Der_nn = (Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i])
//                     + Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i]))/2.;
//
//         double rr_Der_nn = 0.5*(Dx_2_ptpc_2nd(n_v[i+1], n_v[i], n_v[i-1], dr[i]) + Dx_2_ptpc_2nd(n_v_np1[i+1], n_v_np1[i], n_v_np1[i-1], dr[i]) );
//
//         double t_Der_nn = (n_v_np1[i] - n_v[i])/dt;
//
//         double tr_Der_nn = (Dx_ptpc_2nd(n_v_np1[i+1], n_v_np1[i-1], dr[i])
//                             - Dx_ptpc_2nd(n_v[i+1], n_v[i-1], dr[i]))/dt;;
//
//         double r_Der_ss= (Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i])
//                     + Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i]))/2.;
//
//         double rr_Der_ss = 0.5*(Dx_2_ptpc_2nd(s_v[i+1], s_v[i], s_v[i-1], dr[i]) + Dx_2_ptpc_2nd(s_v_np1[i+1], s_v_np1[i], s_v_np1[i-1], dr[i]) );
//
//         double t_Der_ss = (s_v_np1[i] - s_v[i])/dt;
//
//         double tr_Der_ss = (Dx_ptpc_2nd(s_v_np1[i+1], s_v_np1[i-1], dr[i])
//                             - Dx_ptpc_2nd(s_v[i+1], s_v[i-1], dr[i]))/dt;
//
//         gb[i] = get_GB_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//
//         ricci[i] = get_Ricci_Val(r[i], (n_v[i] + n_v_np1[i])/2., r_Der_nn, rr_Der_nn,
//           t_Der_nn, tr_Der_nn,
//           (s_v[i] + s_v_np1[i])/2., r_Der_ss, rr_Der_ss, t_Der_ss, tr_Der_ss);
//       }
//
//       { int i = nx-1;
//
//         gb[i] = 0.;
//         ricci[i] = 0.;
//
//       }
//
//     }
//   }
// }
