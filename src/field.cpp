#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <cmath>
using std::isfinite;
#include <cassert>

#include "grid_data.hpp"
#include "field.hpp"

Field::Field(const std::string name, const std::string type, Grid_data grid):
name{name},
type{type},
v(grid.nx,0.),
grid{grid}
{
}
/*===========================================================================*/
Field::~Field(void)
{
}
/*===========================================================================*/
void Field::set_to_val(const int min, const int max, const double val)
{
   for (int i=min; i<=max; ++i) {
      v[i]= val;
   }
}

/*===========================================================================*/
/*===========================================================================*/
void Field::rescale(){
  assert(v[grid.nx-1] != 0);
  for(int j=0; j< grid.nx; j++){
    v[j] /= v[grid.nx -1];
  }
}
/*===========================================================================*/
void Field::check_field_isfinite(
   const double time, const string level,
   const vector<double> &vec)
{
   size_t index=0; //int index = 0  also works.
   for (auto x: vec) {
      if (!isfinite(x)) {
         cout<<"NaN( ";
         cout<<name<<" ";
         cout<<level<<" ";
         cout<<"time =  "<<time;
         cout<<" ";
         cout<<"index =  "<<index<<" r = "<<grid.r[index];
         cout<<" )"<<endl;
         std::exit(0);
      }
      index++;
   }
}
/*===========================================================================*/
void Field::check_isfinite(const double time)
{
   check_field_isfinite(time,"v",v);
}
/*===========================================================================*/
void Field::check_non_negative(const double time){
  for (int j =0; j<grid.nx; j++){
    if(v[j] + 1e-1< 1e-20){
      cout<<name<<" negative at index = "<<j<<"; value = "<<v[j]<<"; time = "<<time<<endl;
      std::exit(0);
      break;
    }
  }
}
/*===========================================================================*/
/*===========================================================================*/
