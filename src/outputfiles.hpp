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
  void write_anim(const Field &n, const Field &s, const Field &p, const Field &q, const Field &phi);
  void write_fields(const Field &n, const Field &s, const Field &p, const Field &q, const Field &phi);
  void write_characteristics(const std::vector<double> &ingoing, const std::vector<double> &outgoing);
private:
  void write_field(const Field &n);

};

#endif
