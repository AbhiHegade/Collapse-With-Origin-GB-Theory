#ifndef _OUTPUTFILES_HPP_
#define _OUTPUTFILES_HPP_

#include <string>
#include "grid_data.hpp"
#include "field.hpp"

class Write_data{
public:
  Write_data(std::string path);
  ~Write_data(void);

  std::string path;
  void write_grid(const Grid_data grid);
  void write_ah(const Grid_data grid);
  void write_initial_data(const Field &p , const Field &q, const Field &phi);
  void write_fields(const Field &n, const Field &s, const Field &p, const Field &q, const Field &phi);
  void write_characteristics(const std::vector<double> &ingoing, const std::vector<double> &outgoing);
  void write_residual(const std::vector<double> &residual);
  void write_vec(const std::vector<double> &vec, const std::string name);
  void write_MS_mass(const int pos, const double rval,const std::vector<double> &s);
  void write_NER_index(const Grid_data grid);
private:
  void write_field(const Field &n);

};

#endif
