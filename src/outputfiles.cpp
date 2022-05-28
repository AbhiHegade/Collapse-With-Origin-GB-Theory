#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <fstream>

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
  write_output.precision(10);

  for(int i = 0; i<grid.nx-1;i++ ){
    write_output<<grid.x[i]<<",";
  }
  write_output<<grid.x[grid.nx-1]<<"\n";

  write_output.flush();
  write_output.close();

  std::ofstream write_output_1(namer , std::ios::app);
  assert(write_output_1.is_open());
  write_output_1.precision(10);

  for(int i = 0; i<grid.nx-1;i++ ){
    write_output_1<<grid.r[i]<<",";
  }
  write_output_1<<grid.r[grid.nx-1]<<"\n";

  write_output_1.flush();
  write_output_1.close();
}

void Write_data::write_field(const Field &f){
  string name = path + f.name + ".dat";
  std::ofstream write_output(name , std::ios::app);
  assert(write_output.is_open());
  // write_output.setf(std::ios::scientific);
  write_output.precision(10);

  for(int i = 0; i<f.v.size()-1; i++ ){
    write_output<<f.v[i]<<",";
  }
  write_output<<f.v[f.v.size()-1]<<"\n";

  write_output.flush();
  write_output.close();
}

void Write_data::write_fields(const Field &n, const Field &s, const Field &p, const Field &q, const Field &phi){
  write_field(n);
  write_field(s);
  write_field(p);
  write_field(q);
  write_field(phi);

}
