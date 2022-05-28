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
#include "grid_data.hpp"
#include "diagnostics.hpp"
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
  const double ref_val )
  {
    int len = v.size();

    vector<double> v_abs(len);

    for(int i=0; i<len;i++){
      v_abs[i] = fabs(ref_val - v[i]);
      // cout<<"v_abs = "<<v_abs[i]<<endl;
    }

    auto minret = min_element(v_abs.begin(), v_abs.end());

    min_elem = *minret;

    index = minret - v_abs.begin();

    return;
}
//==============================================================================
void Diagnostics::find_apparent_horizon(Grid_data &grid, Field &s_v){

  if(grid.exc_i>0){


  }
  else{
    assert(grid.exc_i==0);
    double min_elem = 0;
    int index = 0;
    const double err_tol= 1e-3;

    find_abs_min(s_v.v,min_elem,index,1.1);
    // cout<<"Minimum = "<<min_elem<<endl;

    if (min_elem<err_tol){
      if (index==0){
        cout<<"Apparent Horizon at the origin."<<endl;
        std::exit(0);
      }
      else{
        cout<<"Found apparent horizon at i = "<<index<<" , "<<"r = "<<grid.r[index]<<" , "<< "t = "<< grid.t_evolve<<endl;
        cout<<"Previous excision point at i = "<<grid.exc_i<<" , "<<"r = "<<grid.r[grid.exc_i]<<endl;
        cout<<"Updating excision point...."<<endl;
        grid.exc_i = index;
        cout<<"done."<<endl;
      }
    }
  }
}
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
