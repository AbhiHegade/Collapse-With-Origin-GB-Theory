#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include <string>
#include <vector>

#include "grid_data.hpp"

/*===========================================================================*/
class Field
{
public:
  const std::string name; // Name of the field
  const std::string type; // Type for KO Dissipation
   std::vector<double> v; //Value of the Field
/*-------------Constructor---------------------------------------------------*/
   Field(const std::string name, const std::string type, Grid_data grid);
/*-------------Destructor----------------------------------------------------*/
   ~Field(void);
/*---------------------------------------------------------------------------*/
   void set_to_val(const int min, const int max, const double val);
   void check_isfinite(const double time);
   void rescale();
   void check_non_negative(const double time);
private:
  Grid_data grid;

   void check_field_isfinite(
      const double time, const std::string message,
      const std::vector<double> &vec);
};

#endif
