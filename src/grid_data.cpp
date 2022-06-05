#include <cmath>
using std::pow;
// using std::fabs;
// #include <iomanip>
// using std::setprecision;
// using std::setw;
// #include <iostream>
// using std::cout;
// using std::endl;
#include <vector>
using std::vector;

#include "grid_data.hpp"
//==============================================================================
Grid_data::Grid_data( const int nx,const int nt, const double l, int exc_i ,const double xl ,const double cl, double t_evolve, const double cfl )
: nx{nx},
  exc_i{exc_i},
  xl{xl},
  cl{cl},
  dx{(cl)/(nx)},
  cfl{cfl},
  x(nx,0.),
  r(nx,0.),
  dr(nx,0.),
  t_evolve{t_evolve},
  dt{cfl*dx},
  nt{nt},
  l{l}
  {
    x[0] = dx/4.;
    r[0] = r_of_x(x[0]);
    dr[0] = r_p_of_x(x[0])*dx;
  for(int j = 1; j<nx-1; j++){
    x[j] = x[0] + (j)*dx;
    r[j] = r_of_x(x[j]);
    dr[j] = r_p_of_x(x[j])*dx;
  }
  x[nx-1] = cl;
  r[nx-1] = 1e100;
  dr[nx-1] = r[nx-1]-r[nx-2];
}
//==============================================================================
Grid_data::~Grid_data(void)
{
}
//==============================================================================

double Grid_data::r_of_x(double x1){

  return (cl*x1)/(cl - x1);
}
//==============================================================================
double Grid_data::r_p_of_x(double x1){
  return pow(cl,2.0)/(pow(cl - x1, 2.0)) ;
}
//==============================================================================
void Grid_data::update_t(){
  t_evolve += dt;
}
//==============================================================================
