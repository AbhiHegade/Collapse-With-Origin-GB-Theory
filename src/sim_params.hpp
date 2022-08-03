#ifndef _READ_PARAM_FILE_HPP_
#define _READ_PARAM_FILE_HPP_

#include<string>

class Sim_params
{
public:
   int nx;
   int nt;
   int save_steps;
   int collapse_and_bh;
   int exc_i;
   double ls;
   double lexp;
   double mu;
   double M;
   double A;
   double rl;
   double ru;
   double cl;


   Sim_params(const std::string output_dir);
   ~Sim_params(void);
private:
   std::string read_sim_params(const std::string output_dir, const std::string find_this_var);
};















#endif
