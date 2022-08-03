#include <cassert>
#include <ctime>
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
using std::to_string;
#include <vector>
using std::vector;
// #include <fstream>
// #include <ctime>
// #include <sys/stat.h>
// #include <sys/types.h>
// #include <sstream>

//==============================================================================

#include "grid_data.hpp"
#include "field.hpp"
#include "solve_metric_fields.hpp"
#include "initial_data.hpp"
#include "diagnostics.hpp"
#include "evolve_scalar_field.hpp"
#include "outputfiles.hpp"
#include "compute_potentials.hpp"
#include "sim_params.hpp"

//==============================================================================

int main(int argc, char const *argv[]) {
  /*-----Creating Directory For Output--------*/

  // std::time_t t_path = std::time(0);
  // std::tm *path_time = std::localtime(&t_path);
  // string path = argv[1];
  // string path = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/" + to_string(path_time->tm_mday) +"_"+to_string(path_time->tm_hour)+ "_"+to_string(path_time->tm_min)+"_"+to_string(path_time->tm_sec) ;
  // char* path_arr;
  // path_arr = &path[0];
  // int rc = mkdir(path_arr, 0755);

  /*----------------------------------------*/
  assert(argc == 2);
  string path = argv[1];

  Sim_params sp(path);
  Write_data write(path);

  //Generate grid with nx, nt, l and exci
  Grid_data grid(sp.nx, sp.nt,sp.ls,sp.lexp,sp.mu,sp.exc_i, sp.cl);
  int save_steps = sp.save_steps;

  //Write grid to file
  write.write_grid(grid);

  //Starting simulation with Minkowski space so M=0
  double M = sp.M;

  //Initialize Evolution classes

  Evolve_scalar_field evolve_scalar_field;
  Solve_metric_fields solve_metric;

  //Initialize Initial_data class
  Initial_data initialdata(sp.A,sp.ru,sp.rl, M);

  //Initialize diagonistics class
  Diagnostics diagnostics;

  //============================================================================
  //Initializing Field objects for the simulation

  Field s("shift","odd",grid);
  Field n("lapse", "even", grid);
  Field p("dt_phi", "even", grid);
  Field q("dr_phi", "odd", grid);
  Field phi("phi", "even", grid);


  //============================================================================
  vector<double> residual(grid.nx,0 );
  vector<double> gb(grid.nx,0);
  vector<double> n_nm1(grid.nx,0);
  vector<double> s_nm1(grid.nx,0);
  vector<double> p_nm1(grid.nx,0);
  vector<double> q_nm1(grid.nx,0);
  vector<double> phi_nm1(grid.nx,0);
  vector<double> ingoing(grid.nx,0);
  vector<double> outgoing(grid.nx,0);

  //============================================================================
  /* Simulation Parameters */

  // cout<<"Simulation Parameters"<<"\n"<<endl;
  // cout<<"---------------------------------------------------------------"<<endl;
  // cout<<"---------------------------------------------------------------"<<endl;
  cout<<"Saving file at : "<<path<<endl;
  cout<<"nx = "<<grid.nx<<endl;
  cout<<"nt = "<<grid.nt<<endl;
  cout<<"cl = "<<grid.cl<<endl;
  cout<<"t_save_steps = "<<save_steps<<endl;
  cout<<"dx = "<<grid.dx<<endl;
  cout<<"dt = "<<grid.dt<<endl;
  cout<<"shiftsymm coupling ls = "<<grid.ls<<endl;
  cout<<"exp coupling lexp = "<<grid.lexp<<endl;
  cout<<"mu coupling mu = "<<grid.mu<<endl;
  cout<<"Scalar Field Amplitude = "<<initialdata.amp<<endl;
  cout<<"ru = "<<initialdata.r_u<<endl;
  cout<<"rl = "<<initialdata.r_l<<endl;
  cout<<"Starting simulation..."<<endl;

  //============================================================================
  //Set initial data
  if(M<fabs(1e-16)){
  initialdata.set_Minkowski(grid, n, s, p, q, phi);
  }
  else{
    initialdata.set_bh_bump(grid,n,s,p,q,phi );
  }
  //Check if initial data is a naked singularity
  {
  vector<double> ns_check(grid.nx,0);
  for(int j =grid.exc_i; j<grid.nx; j++){
    ns_check[j] = grid.r[j] - 8.*q.v[j]*beta_p(grid.ls,grid.lexp,grid.mu, phi.v[j]);
  }
  {
    double min_elem = 0.;
    int start_index = 100;
    int index = 0;

  diagnostics.find_abs_min(ns_check, min_elem, index, 0, start_index);
  if((min_elem<1e-2) && (index > 100)){
    cout<<"Naked sing check val = "<<min_elem<<endl;
    cout<<"Data too strong leads to naked singularity."<<endl;
    cout<<"NaN at index = "<<index<<endl;
  }
  }
  }
  //===========================================================================
  write.write_initial_data(p,q,phi);
  //===========================================================================
  //===========================================================================
  //Check if apparent horizon is present in initial data (only relevant for black hole initial data).
  diagnostics.find_apparent_horizon(grid,s);

  //Solve for metric fields
  solve_metric.solve( grid, n, s , p ,q,phi);

  int mass_extraction_radius = 3*(grid.nx/4);
  cout<<"Initial MS_mass = "<< setprecision(4)<<grid.r[mass_extraction_radius]*(pow((s.v[mass_extraction_radius]),2.)/2.)<<endl;
  // std::exit(0);
  //==========================================
  diagnostics.find_apparent_horizon(grid,s);
  cout<< "collapse_and_bh = "<< sp.collapse_and_bh<< endl;
  if(sp.collapse_and_bh == 0){
    if(grid.exc_i>0){
      cout<<"BH formation at t=0."<<endl;
      std::exit(0);
    }
    else{
      cout<<"No BH formation at t=0."<<endl;
      std::exit(0);
    }
  }
  // diagnostics.check_for_elliptic_region(grid, n, s, p, q, phi, ingoing, outgoing);

  //Write data to file
  write.write_fields(n, s , p ,q, phi);
  write.write_characteristics(ingoing, outgoing);


  int i_e = 0;
  while(i_e<grid.nt){

  n_nm1 = n.v;
  s_nm1 = s.v;
  p_nm1 = p.v;
  q_nm1 = q.v;
  phi_nm1 = phi.v;

  evolve_scalar_field.evolve(grid, n, s, p, q, phi);
  grid.update_t();
  solve_metric.solve( grid, n, s , p ,q,phi);
  diagnostics.compute_e_rr_residual(grid, n_nm1, s_nm1, p_nm1, q_nm1, phi_nm1, s.v, p.v,residual);
  diagnostics.find_apparent_horizon(grid,s);
  diagnostics.check_for_elliptic_region(grid, n, s, p, q, phi, ingoing, outgoing);

  i_e += 1;

  if ((i_e%save_steps ==0) ){
    diagnostics.compute_GB(grid,
    n_nm1,
    s_nm1,
    s.v,
    n.v,
    gb);
    write.write_vec(gb, "gb");
    write.write_residual(residual);
    write.write_fields(n, s , p ,q, phi);
    write.write_characteristics(ingoing, outgoing);
    write.write_ah(grid);

  }


}

  cout<<"Final time = "<<grid.t_evolve<<endl;
  if(grid.exc_i>0){
    int mass_extraction_radius = 3*(grid.nx/4);
    cout<<"exit_code_1, BH_Formation, MS_mass = "<< setprecision(4)<<grid.r[mass_extraction_radius]*(pow((s.v[mass_extraction_radius]),2.)/2.)<<", run finished successfully."<<endl;
  }
  else{
    cout<<"exit_code_0, no black hole formation, MS_mass = "<< setprecision(4)<<grid.r[mass_extraction_radius]*(pow((s.v[mass_extraction_radius]),2.)/2.)<<", run finished successfully."<<endl;
  }
  // cout<<"---------------------------------------------------------------"<<endl;
  // cout<<"---------------------------------------------------------------"<<endl;
  return 0;
}
