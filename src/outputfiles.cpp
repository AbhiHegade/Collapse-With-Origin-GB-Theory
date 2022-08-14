#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <fstream>
#include <vector>
using std::vector;
#include<cassert>

#include "outputfiles.hpp"
#include "grid_data.hpp"
#include "field.hpp"

Write_data::Write_data(string path):
path{path+"/"}{

}

Write_data::~Write_data(void){

}

void Write_data::write_grid(const Grid_data grid){
  string namex = path + "x" + ".dat";
  string namer = path + "r" + ".dat";
  std::ofstream write_output(namex , std::ios::app);
  assert(write_output.is_open());
  write_output.precision(16);

  for(int i = 0; i<grid.nx-1;i++ ){
    write_output<<grid.x[i]<<",";
  }
  write_output<<grid.x[grid.nx-1]<<"\n";

  write_output.flush();
  write_output.close();

  std::ofstream write_output_1(namer , std::ios::app);
  assert(write_output_1.is_open());
  write_output_1.precision(16);

  for(int i = 0; i<grid.nx-1;i++ ){
    write_output_1<<grid.r[i]<<",";
  }
  write_output_1<<grid.r[grid.nx-1]<<"\n";

  write_output_1.flush();
  write_output_1.close();
}

void Write_data::write_ah(const Grid_data grid){
  string name_ah = path + "ah" + ".dat";
  string name_exci = path + "exci.dat";
  std::ofstream write_output(name_ah , std::ios::app);
  assert(write_output.is_open());
  write_output.precision(16);
  write_output<<grid.ah_index<<" ";
  write_output.flush();
  write_output.close();

  std::ofstream write_output_1(name_exci , std::ios::app);
  assert(write_output_1.is_open());
  write_output_1.precision(16);
  write_output_1<<grid.exc_i<<" ";
  write_output_1.flush();
  write_output_1.close();
}
void Write_data::write_field(const Field &f){
  string name = path + f.name + ".dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  // write_output.setf(std::ios::scientific);
  write_output.precision(16);

  for(int i = 0; i<f.v.size()-1; i++ ){
    write_output<<f.v[i]<<",";
  }
  write_output<<f.v[f.v.size()-1]<<"\n";

  write_output.flush();
  write_output.close();
}
void Write_data::write_initial_data(const Field &p, const Field &q,const Field &phi){
  string name = path + "p_init"+ ".dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  // write_output.setf(std::ios::scientific);
  write_output.precision(16);

  for(int i = 0; i<p.v.size()-1; i++ ){
    write_output<<p.v[i]<<",";
  }
  write_output<<p.v[p.v.size()-1]<<"\n";

  write_output.flush();
  write_output.close();

  string name_q = path + "q_init"+ ".dat";
  std::ofstream write_output_1(name_q , std::ios::app);
  assert(write_output_1.is_open());
  // write_output.setf(std::ios::scientific);
  write_output_1.precision(16);

  for(int i = 0; i<q.v.size()-1; i++ ){
    write_output_1<<q.v[i]<<",";
  }
  write_output_1<<q.v[q.v.size()-1]<<"\n";

  write_output_1.flush();
  write_output_1.close();

  string name_phi = path + "phi_init"+ ".dat";
  std::ofstream write_output_2(name_phi , std::ios::app);
  assert(write_output_2.is_open());
  // write_output.setf(std::ios::scientific);
  write_output_2.precision(16);

  for(int i = 0; i<phi.v.size()-1; i++ ){
    write_output_2<<phi.v[i]<<",";
  }
  write_output_2<<phi.v[phi.v.size()-1]<<"\n";

  write_output_2.flush();
  write_output_2.close();
}
void Write_data::write_fields(const Field &n, const Field &s, const Field &p, const Field &q, const Field &phi){
  write_field(n);
  write_field(s);
  write_field(p);
  write_field(q);
  write_field(phi);

}

void Write_data::write_characteristics(const vector<double> &ingoing, const vector<double> &outgoing){
  write_vec(ingoing, "ingoing");
  write_vec(outgoing, "outgoing");


}

void Write_data::write_residual(const vector<double> &residual){
  write_vec(residual, "e_rr_res");

}

void Write_data::write_vec(const vector<double> &vec, const string name_v){
  string name = path + name_v + ".dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  write_output.precision(16);

  int nx = vec.size();
  for(int i = 0; i<nx-1; i++ ){
    write_output<<vec[i]<<",";
  }
  write_output<<vec[nx-1]<<"\n";

  write_output.flush();
  write_output.close();

}

void Write_data::write_MS_mass(const int pos, const double rval, const vector<double> &s){
  string name = path+"ms_mass.dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  write_output.precision(16);

  double MS_mass = 0.5*rval*(s[pos]*s[pos]);
  write_output<<MS_mass<<"\n";
  write_output.flush();
  write_output.close();
}
