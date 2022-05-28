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

#include "field.hpp"
#include "initial_data.hpp"
#include "grid_data.hpp"
//==============================================================================
Initial_data::Initial_data(const double amp, const double ru, const double rl,const double M):
M{M},
amp{amp},
r_u{ru},
r_l{rl}
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
        p_v.v[i]= 0;
        q_v.v[i]= (
        2*(r-r_l)*pow(r_u-r,2)
        -	2*pow(r-r_l,2)*(r_u-r)
        +	pow(r_u-r,2)
        -	pow(r-r_l,2)
        )*bump;
     } else {
        phi.v[i]= 0;
        p_v.v[i]= 0;
        q_v.v[i]= 0;
     }
     max_vN = (fabs(phi.v[i])>max_vN) ? fabs(phi.v[i]) : max_vN;

     n_v.v[i]= 1.0;
     s_v.v[i]= 1e-20;
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
// void Initial_data::set_Minkowski(Grid_data grid, Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi_v){
//   assert(fabs(M)<1e-16);
//   assert(r_u>r_l);
//   assert(grid.exc_i==0);
//   double max_vN= 0;
//   for (int i=grid.exc_i; i<grid.nx-1; ++i) {
//      double r= grid.r[i];
//      phi_v.v[i] = 0.01*exp(-pow((r - 8.),2.)/(2.*1.));
//      q_v.v[i] = 0.01*exp(-pow((r - 8.),2.)/(2.*1.))*(-(r-8.));
//      p_v.v[i] = 0.;
//      n_v.v[i] = 1.;
//      s_v.v[i] = 1e-20;
//   }
//   {
//     int i = grid.nx-1;
//     phi_v.v[i] = 0.;
//     q_v.v[i] = 0.;
//     p_v.v[i] = 0.;
//     n_v.v[i] = 1.;
//     s_v.v[i] = 1e-20;
//
//   }
//   /* rescale so amp is actual maximum val */
// }
//==============================================================================
void Initial_data::set_bh_bump(Grid_data grid,  Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi){
  cout<<"Not implemented"<<endl;
  std::exit(0);
//   assert(M>0);
//   assert(r_u>r_l);
//   assert(grid.exc_i>0);
//   double max_vN= 0;
//   for (int i=grid.exc_i; i<grid.nx-1; ++i) {
//      double r= grid.r[i];
//      if ((r>r_l)
//      &&  (r<r_u)
//      ) {
//         double bump= exp(-1./(r_u-r))*exp(-1./(r-r_l));
//
//         phi.v[i]= pow(r-r_l,2)*pow(r_u-r,2)*bump;
//         p_v.v[i]= 0.0;
//         q_v.v[i]= (
//         2*(r-r_l)*pow(r_u-r,2)
//         -	2*pow(r-r_l,2)*(r_u-r)
//         +	pow(r_u-r,2)
//         -	pow(r-r_l,2)
//         )*bump;
//      } else {
//         phi.v[i]= 0;
//         p_v.v[i]= 0;
//         q_v.v[i]= 0;
//      }
//      max_vN = (fabs(phi.v[i])>max_vN) ? fabs(phi.v[i]) : max_vN;
//
//      n_v.v[i]= 1.0;
//      s_v.v[i]= pow(2*M/r,0.5);
//   }
//   /* rescale so amp is actual maximum val */
//   for (int i=0; i<grid.nx-1; ++i) {
//      phi.v[i]*= amp/max_vN;
//      q_v.v[i]*= amp/max_vN;
//      p_v.v[i]*= amp/max_vN;
//   }
//   /* left boundary values */
//   s_v.lbval = pow(2*M/grid.r[grid.exc_i],0.5);
//   /* right boundary values */
//   phi.rbval = 0;
//   p_v.rbval = 0;
//   q_v.rbval = 0;
//   n_v.rbval = 1.0;
//   s_v.rbval = 0;
// /* Update Derivatives*/
// s_v.compute_r_Der_ord_4();
// n_v.compute_r_Der_ord_4();
// p_v.compute_r_Der_ord_4();
// q_v.compute_r_Der_ord_4();
// phi.compute_r_Der_ord_4();
}

//==============================================================================
