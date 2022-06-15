#include <cassert>
#include <ctime>
#include <cmath>
#include <string>
using std::string;
#include <iomanip>
using std::setprecision;
using std::setw;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
using std::to_string;
#include <vector>
using std::vector;
#include <fstream>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <sstream>

//==============================================================================
#include "compute_potentials.hpp"
#include "grid_data.hpp"
#include "field.hpp"
#include "solve_metric_fields.hpp"
#include "initial_data.hpp"
#include "diagnostics.hpp"
#include "evolve_scalar_field.hpp"
#include "outputfiles.hpp"


int main(int argc, char const *argv[]) {
  /*-----Creating Directory For Output--------*/
  std::time_t t_path = std::time(0);
  std::tm *path_time = std::localtime(&t_path);

  string path = "/Users/abhi/Work/Projects/Hyperbolitcity-Gravitational-Collapse/code-f-phi/output/Convergence-Runs/" + to_string(path_time->tm_mday) +"_"+to_string(path_time->tm_hour)+ "_"+to_string(path_time->tm_min)+"_"+to_string(path_time->tm_sec);
  string path_h = path + "_h";
  string path_2h = path + "_2h";
  string path_4h = path + "_4h";
  char* path_arr_h;
  char* path_arr_2h;
  char* path_arr_4h;

  path_arr_h = &path_h[0];
  path_arr_2h = &path_2h[0];
  path_arr_4h = &path_4h[0];
  int rc_h = mkdir(path_arr_h, 0755), rc_2h = mkdir(path_arr_2h, 0755), rc_4h = mkdir(path_arr_4h, 0755);
  /*----------------------------------------*/
  Write_data write_h(path_h);
  Write_data write_2h(path_2h);
  Write_data write_4h(path_4h);

  int nx_h = 8000, nt_h = 4000;
  double l = 0.01;
  int save_steps_h = 20;
  int update_ah_h = 8;

  assert((nx_h%4 == 0)&&(nt_h%4 == 0));
  assert(save_steps_h%4 == 0);
  assert(update_ah_h%4==0);

  Grid_data grid_h(nx_h, nt_h,l);
  Grid_data grid_2h(nx_h/2, nt_h/2, l);
  Grid_data grid_4h(nx_h/4, nt_h/4, l);


  int save_steps_2h = save_steps_h/2;
  int save_steps_4h = save_steps_h/4;

  int update_ah_2h = update_ah_h/2;
  int update_ah_4h = update_ah_h/4;

  write_h.write_grid(grid_h);
  write_2h.write_grid(grid_2h);
  write_4h.write_grid(grid_4h);

  //Initialize Evolution classes
  Evolve_scalar_field evolve_scalar_field;
  Solve_metric_fields solve_metric;
  //Initialize Initial_data class
  double M = 0.;
  Initial_data initialdata(0.27,12.,8., M);
  Diagnostics diagnostics;

  //============================================================================

  //============================================================================
  /* Simulation Parameters */

  cout<<"Simulation Parameters"<<"\n"<<"---------------------------------------"<<endl;
  cout<<"Run type = Convergence checks"<<endl;
  cout<<"Saving file at : path = "<<"\'"+path+"\'"<<endl;
  cout<<"nx_h = "<<grid_h.nx<<endl;
  cout<<"nt_h = "<<grid_h.nt<<endl;
  cout<<"t_save_steps_h = "<<save_steps_h<<endl;
  cout<<"dx = "<<grid_h.dx<<endl;
  cout<<"dt = "<<grid_h.dt<<endl;
  cout<<"GB coupling l = "<<grid_h.l<<endl;
  cout<<"Scalar Field Amplitude = "<<initialdata.amp<<endl;
  cout<<"Scalar Field ru = "<<initialdata.r_u<<" , rl = "<<initialdata.r_l<<endl;
  cout<<"Starting simulation..."<<endl;

  //============================================================================


{

  //Initializing Field objects for the simulation
  Field s_4h("shift","odd",grid_4h);
  Field n_4h("lapse", "even", grid_4h);
  Field p_4h("dt_phi", "even", grid_4h);
  Field q_4h("dr_phi", "odd", grid_4h);
  Field phi_4h("phi", "even", grid_4h);

  vector<double> residual_4h(grid_4h.nx,0 );
  vector<double> n_nm1_4h(grid_4h.nx,0);
  vector<double> s_nm1_4h(grid_4h.nx,0);
  vector<double> p_nm1_4h(grid_4h.nx,0);
  vector<double> q_nm1_4h(grid_4h.nx,0);
  vector<double> phi_nm1_4h(grid_4h.nx,0);
  //============================================================================
  vector<double> ingoing_4h(grid_4h.nx,0);
  vector<double> outgoing_4h(grid_4h.nx,0);

  initialdata.set_Minkowski(grid_4h,  n_4h, s_4h, p_4h, q_4h, phi_4h); //Set initial data

  diagnostics.find_apparent_horizon(grid_4h,s_4h); //Check if apparent horizon is present in initial data
  solve_metric.solve( grid_4h, n_4h, s_4h , p_4h ,q_4h ,phi_4h); //Solve for metric fields
  write_4h.write_fields(n_4h, s_4h , p_4h ,q_4h , phi_4h); //Write data to file
  diagnostics.find_apparent_horizon(grid_4h,s_4h); //Check if apparent horizon is present
  diagnostics.check_for_elliptic_region(grid_4h, n_4h, s_4h, p_4h, q_4h, phi_4h, ingoing_4h, outgoing_4h);
  write_4h.write_characteristics(ingoing_4h, outgoing_4h);

  int i_e = 0;
  while(i_e<grid_4h.nt){

  n_nm1_4h = n_4h.v;
  s_nm1_4h = s_4h.v;
  p_nm1_4h = p_4h.v;
  q_nm1_4h = q_4h.v;
  phi_nm1_4h = phi_4h.v;

  evolve_scalar_field.evolve(grid_4h, n_4h, s_4h, p_4h, q_4h, phi_4h);


  grid_4h.update_t();

  solve_metric.solve( grid_4h, n_4h, s_4h , p_4h ,q_4h, phi_4h);
  diagnostics.compute_e_rr_residual(grid_4h, n_nm1_4h, s_nm1_4h, p_nm1_4h, q_nm1_4h, phi_nm1_4h, s_4h.v, p_4h.v,residual_4h);
  diagnostics.find_apparent_horizon(grid_4h,s_4h);
  diagnostics.check_for_elliptic_region(grid_4h, n_4h, s_4h, p_4h, q_4h, phi_4h, ingoing_4h, outgoing_4h);
  i_e += 1;
  if ((i_e%save_steps_4h ==0) ){
    // cout<<"Saving file now nt = "<<i_e<<endl;
    write_4h.write_residual(residual_4h);
    write_4h.write_fields(n_4h, s_4h , p_4h ,q_4h, phi_4h);
    write_4h.write_characteristics(ingoing_4h, outgoing_4h);
  }

}


}
  cout<<"Level 4h done."<<endl;

{
    //Initializing Field objects for the simulation
    Field s_2h("shift","odd",grid_2h);
    Field n_2h("lapse", "even", grid_2h);
    Field p_2h("dt_phi", "even", grid_2h);
    Field q_2h("dr_phi", "odd", grid_2h);
    Field phi_2h("phi", "even", grid_2h);

    vector<double> residual_2h(grid_2h.nx,0 );
    vector<double> n_nm1_2h(grid_2h.nx,0);
    vector<double> s_nm1_2h(grid_2h.nx,0);
    vector<double> p_nm1_2h(grid_2h.nx,0);
    vector<double> q_nm1_2h(grid_2h.nx,0);
    vector<double> phi_nm1_2h(grid_2h.nx,0);
    //============================================================================
    vector<double> ingoing_2h(grid_2h.nx,0);
    vector<double> outgoing_2h(grid_2h.nx,0);

    initialdata.set_Minkowski(grid_2h,  n_2h, s_2h, p_2h, q_2h, phi_2h); //Set initial data

    diagnostics.find_apparent_horizon(grid_2h,s_2h); //Check if apparent horizon is present in initial data
    solve_metric.solve( grid_2h, n_2h, s_2h , p_2h ,q_2h ,phi_2h); //Solve for metric fields
    write_2h.write_fields(n_2h, s_2h , p_2h ,q_2h , phi_2h); //Write data to file
    diagnostics.find_apparent_horizon(grid_2h,s_2h); //Check if apparent horizon is present
    diagnostics.check_for_elliptic_region(grid_2h, n_2h, s_2h, p_2h, q_2h, phi_2h, ingoing_2h, outgoing_2h);
    write_2h.write_characteristics(ingoing_2h, outgoing_2h);

    int i_e = 0;
    while(i_e<grid_2h.nt){

    n_nm1_2h = n_2h.v;
    s_nm1_2h = s_2h.v;
    p_nm1_2h = p_2h.v;
    q_nm1_2h = q_2h.v;
    phi_nm1_2h = phi_2h.v;

    evolve_scalar_field.evolve(grid_2h, n_2h, s_2h, p_2h, q_2h, phi_2h);


    grid_2h.update_t();

    solve_metric.solve( grid_2h, n_2h, s_2h , p_2h ,q_2h, phi_2h);
    diagnostics.compute_e_rr_residual(grid_2h, n_nm1_2h, s_nm1_2h, p_nm1_2h, q_nm1_2h, phi_nm1_2h, s_2h.v, p_2h.v,residual_2h);
    diagnostics.find_apparent_horizon(grid_2h,s_2h);
    diagnostics.check_for_elliptic_region(grid_2h, n_2h, s_2h, p_2h, q_2h, phi_2h, ingoing_2h, outgoing_2h);
    i_e += 1;
    if ((i_e%save_steps_2h ==0) ){
      // cout<<"Saving file now nt = "<<i_e<<endl;
      write_2h.write_residual(residual_2h);
      write_2h.write_fields(n_2h, s_2h , p_2h ,q_2h, phi_2h);
      write_2h.write_characteristics(ingoing_2h, outgoing_2h);
    }

  }


}
  cout<<"Level 2h done."<<endl;

{
  //Initializing Field objects for the simulation
  Field s_h("shift","odd",grid_h);
  Field n_h("lapse", "even", grid_h);
  Field p_h("dt_phi", "even", grid_h);
  Field q_h("dr_phi", "odd", grid_h);
  Field phi_h("phi", "even", grid_h);

  vector<double> residual_h(grid_h.nx,0 );
  vector<double> n_nm1_h(grid_h.nx,0);
  vector<double> s_nm1_h(grid_h.nx,0);
  vector<double> p_nm1_h(grid_h.nx,0);
  vector<double> q_nm1_h(grid_h.nx,0);
  vector<double> phi_nm1_h(grid_h.nx,0);
  //============================================================================
  vector<double> ingoing_h(grid_h.nx,0);
  vector<double> outgoing_h(grid_h.nx,0);

  initialdata.set_Minkowski(grid_h,  n_h, s_h, p_h, q_h, phi_h); //Set initial data

  diagnostics.find_apparent_horizon(grid_h,s_h); //Check if apparent horizon is present in initial data
  solve_metric.solve( grid_h, n_h, s_h , p_h ,q_h ,phi_h); //Solve for metric fields
  write_h.write_fields(n_h, s_h , p_h ,q_h , phi_h); //Write data to file
  diagnostics.find_apparent_horizon(grid_h,s_h); //Check if apparent horizon is present
  diagnostics.check_for_elliptic_region(grid_h, n_h, s_h, p_h, q_h, phi_h, ingoing_h, outgoing_h);
  write_h.write_characteristics(ingoing_h, outgoing_h);

  int i_e = 0;
  while(i_e<grid_h.nt){

  n_nm1_h = n_h.v;
  s_nm1_h = s_h.v;
  p_nm1_h = p_h.v;
  q_nm1_h = q_h.v;
  phi_nm1_h = phi_h.v;

  evolve_scalar_field.evolve(grid_h, n_h, s_h, p_h, q_h, phi_h);


  grid_h.update_t();

  solve_metric.solve( grid_h, n_h, s_h , p_h ,q_h, phi_h);
  diagnostics.compute_e_rr_residual(grid_h, n_nm1_h, s_nm1_h, p_nm1_h, q_nm1_h, phi_nm1_h, s_h.v, p_h.v,residual_h);
  diagnostics.find_apparent_horizon(grid_h,s_h);
  diagnostics.check_for_elliptic_region(grid_h, n_h, s_h, p_h, q_h, phi_h, ingoing_h, outgoing_h);
  i_e += 1;
  if ((i_e%save_steps_h ==0) ){
    // cout<<"Saving file now nt = "<<i_e<<endl;
    write_h.write_residual(residual_h);
    write_h.write_fields(n_h, s_h , p_h ,q_h, phi_h);
    write_h.write_characteristics(ingoing_h, outgoing_h);
  }
}


}
cout<<"Level h done."<<endl;
  cout<<"Run finished successfully"<<endl;
  cout<<"Final time = "<<grid_h.t_evolve<<endl;
  return 0;
}
