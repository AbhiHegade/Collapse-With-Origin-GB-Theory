#ifndef _GRID_DATA_HPP_
#define _GRID_DATA_HPP_

#include <cmath>
#include <vector>

class Grid_data{
public:
  /*----Constructor-----------------------------------------------------------*/
  Grid_data(const int nx,const int nt, const double ls, const double lexp, const double mu,
    int exc_i = 0 ,const double cl = 100.0,const double xl = 0., double t_evolve = 0.0, const double cfl = 0.2 );
  /*----Destructor------------------------------------------------------------*/
  ~Grid_data(void);
  //----------------------------------------------------------------------------
  const int nx; //Number of grid points in x
  int exc_i; //Excision location
  const double xl; //Left end of x grid
  const double cl; // Compactification Length
  const double dx; // dx
  const double cfl;
  int ah_index;

  std::vector<double> x; // x grid values stored in a vector
  std::vector<double> r; // r grid values stored in a vector
  std::vector<double> dr;// dr(x) stored in a vector
  double t_evolve; // Evolution Time
  const double dt;
  const int nt;
  const double ls; //Coupling constant for shift-symmetric term
  const double lexp; //Coupling for exponential term
  const double mu; // Strength of exponential coupling

  void update_t();

private:

};

#endif
