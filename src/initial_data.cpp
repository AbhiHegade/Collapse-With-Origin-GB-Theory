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
#include <cassert>

#include "field.hpp"
#include "initial_data.hpp"
#include "grid_data.hpp"
//==============================================================================
Initial_data::Initial_data(const double amp, const double r_u, const double r_l,const double M):
M{M},
amp{amp},
r_u{r_u},
r_l{r_l}
{

}
//==============================================================================
Initial_data::~Initial_data(void){

}
//==============================================================================
void Initial_data::set_Minkowski(Grid_data grid, Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi){
  assert(fabs(M)<1e-16);
  assert(r_u>r_l);
  assert(grid.exc_i==0);
  double max_vN= 0;
  for (int i=grid.exc_i; i<grid.nx-1; ++i) {
     double r= grid.r[i];
     if ((r>r_l)
     &&  (r<r_u)
     ) {
        double bump= exp(-1./(r_u-r))*exp(-1./(r-r_l));

        phi.v[i]= pow(r-r_l,2)*pow(r_u-r,2)*bump;
        q_v.v[i]= (
        2*(r-r_l)*pow(r_u-r,2)
        -	2*pow(r-r_l,2)*(r_u-r)
        +	pow(r_u-r,2)
        -	pow(r-r_l,2)
        )*bump;
        p_v.v[i]= q_v.v[i] + (phi.v[i]/r);
        // p_v.v[i] = 0.;
     } else {
        phi.v[i]= 0;
        p_v.v[i]= 0;
        q_v.v[i]= 0;
     }
     max_vN = (fabs(phi.v[i])>max_vN) ? fabs(phi.v[i]) : max_vN;

     n_v.v[i]= 1.0;
     s_v.v[i]= 1e-10;
  }
  /* rescale so amp is actual maximum val */
  for (int i=0; i<grid.nx-1; ++i) {
     phi.v[i]*= amp/max_vN;
     q_v.v[i]*= amp/max_vN;
     p_v.v[i]*= amp/max_vN;
  }
  {
    int i = grid.nx -1;
    n_v.v[i] = 1.0;
    s_v.v[i] = 0.;
    phi.v[i] = 0.;
    q_v.v[i] = 0.;
    p_v.v[i] = 0.;

  }
}
//==============================================================================
void Initial_data::set_bh_bump(Grid_data grid,  Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi){
  // cout<<"Not implemented"<<endl;
  // std::exit(0);
  assert(M>0);
  assert(r_u>r_l);
  assert(grid.exc_i>0);
  assert(grid.r[grid.exc_i]<2.*M);
  double max_vN= 0;
  for (int i=grid.exc_i; i<grid.nx-1; ++i) {
     double r= grid.r[i];
     if ((r>r_l)
     &&  (r<r_u)
     ) {
        double bump= exp(-1./(r_u-r))*exp(-1./(r-r_l));

        phi.v[i]= pow(r-r_l,2)*pow(r_u-r,2)*bump;
        q_v.v[i]= (
        2*(r-r_l)*pow(r_u-r,2)
        -	2*pow(r-r_l,2)*(r_u-r)
        +	pow(r_u-r,2)
        -	pow(r-r_l,2)
        )*bump;
        p_v.v[i]= q_v.v[i] + (phi.v[i]/r);
        // p_v.v[i]= 0.0;
     } else {
        phi.v[i]= 0;
        p_v.v[i]= 0;
        q_v.v[i]= 0;
     }
     max_vN = (fabs(phi.v[i])>max_vN) ? fabs(phi.v[i]) : max_vN;

     n_v.v[i]= 1.;
     s_v.v[i]= pow(2.*M/r,0.5);
  }
  /* rescale so amp is actual maximum val */
  for (int i=0; i<grid.nx-1; ++i) {
     phi.v[i]*= amp/max_vN;
     q_v.v[i]*= amp/max_vN;
     p_v.v[i]*= amp/max_vN;
  }
  {
    int i = grid.nx -1;
    n_v.v[i] = 1.0;
    s_v.v[i] = 0.;
    phi.v[i] = 0.;
    q_v.v[i] = 0.;
    p_v.v[i] = 0.;

  }
}

//==============================================================================
void Initial_data::set_Gaussian(Grid_data grid,  Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi){
  // cout<<"Not implemented"<<endl;
  // std::exit(0);
  if(M>1e-20){
  assert(r_u>r_l);
  assert(grid.exc_i>0);
  assert(grid.r[grid.exc_i]<2.*M);
  double w0 = (r_u - r_l)/2.;
  double r0 = (r_l+r_u)/2.;
  double max_vN= 0;
  for (int i=grid.exc_i; i<grid.nx-1; ++i) {
    double r = grid.r[i];
    double exparg = (r-r0)/w0;
    double bump = exp(- (exparg*exparg));
    phi.v[i] = pow((r/w0),2.)*bump;
    q_v.v[i] = (2.*bump*r/pow(w0,4.))*(pow(w0,2.) - pow(r,2.) + r*r0);
    p_v.v[i]= q_v.v[i] + (phi.v[i]/r);
    max_vN = (fabs(phi.v[i])>max_vN) ? fabs(phi.v[i]) : max_vN;

   n_v.v[i]= 1.;
   s_v.v[i]= pow(2.*M/r,0.5);
  }
  /* rescale so amp is actual maximum val */
  for (int i=0; i<grid.nx-1; ++i) {
     phi.v[i]*= amp/max_vN;
     q_v.v[i]*= amp/max_vN;
     p_v.v[i]*= amp/max_vN;
  }
  {
    int i = grid.nx -1;
    n_v.v[i] = 1.0;
    s_v.v[i] = 0.;
    phi.v[i] = 0.;
    q_v.v[i] = 0.;
    p_v.v[i] = 0.;

  }
}
else{
  assert(r_u>r_l);
  assert(grid.exc_i==0);
  double w0 = (r_u - r_l)/2.;
  double r0 = (r_l+r_u)/2.;
  double max_vN= 0;
  for (int i=grid.exc_i; i<grid.nx-1; ++i) {
    double r = grid.r[i];
    double exparg = (r-r0)/w0;
    double bump = exp(- (exparg*exparg));
    phi.v[i] = pow((r/w0),2.)*bump;
    q_v.v[i] = (2.*bump*r/pow(w0,4.))*(pow(w0,2.) - pow(r,2.) + r*r0);
    p_v.v[i]= q_v.v[i] + (phi.v[i]/r);

    max_vN = (fabs(phi.v[i])>max_vN) ? fabs(phi.v[i]) : max_vN;

   n_v.v[i]= 1.;
   s_v.v[i]= 1e-8;
  }
  /* rescale so amp is actual maximum val */
  for (int i=0; i<grid.nx-1; ++i) {
     phi.v[i]*= amp/max_vN;
     q_v.v[i]*= amp/max_vN;
     p_v.v[i]*= amp/max_vN;
  }
  {
    int i = grid.nx -1;
    n_v.v[i] = 1.0;
    s_v.v[i] = 0.;
    phi.v[i] = 0.;
    q_v.v[i] = 0.;
    p_v.v[i] = 0.;

  }

}
}
