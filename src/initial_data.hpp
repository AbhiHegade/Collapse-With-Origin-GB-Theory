#ifndef _INITIAL_DATA_HPP_
#define _INITIAL_DATA_HPP_

#include <cmath>
#include <vector>
#include <string>


#include "field.hpp"
#include "grid_data.hpp"

class Initial_data{
public:
  Initial_data(const double amp, const double r_u, const double r_l,const double M = 0);
  ~Initial_data(void);
  const double M;
  const double amp;
  const double r_u;
  const double r_l;
  void set_Minkowski(Grid_data grid,Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi);
  void set_bh_bump(Grid_data grid, Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi);
  void set_Gaussian(Grid_data grid, Field &n_v, Field &s_v, Field &p_v, Field &q_v, Field &phi);
};


#endif
