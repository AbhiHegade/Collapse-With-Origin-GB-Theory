#include <cassert>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <string>
using std::string;
#include <iomanip>
#include <fstream>

#include "sim_params.hpp"
/*===========================================================================*/
string Sim_params::read_sim_params(const string output_dir, const string find_this_var)
{
   ifstream infile(output_dir+"/sim_params.txt");
   assert(infile.good());

   string name, val;

   while (infile >> name >> val) {
      if (name==find_this_var) {
         return val;
      }
   }
   cout << "ERROR: did not find " << find_this_var << endl;
   std::exit(0);
}
/*===========================================================================*/
Sim_params::Sim_params(const string output_dir)
{
   nx = stoi(read_sim_params(output_dir,"nx"));
   nt= stoi(read_sim_params(output_dir,"nt"));
   save_steps= stoi(read_sim_params(output_dir,"save_steps"));
   l= stod(read_sim_params(output_dir,"l"));
   M = stod(read_sim_params(output_dir, "initial_mass"));
   A = stod(read_sim_params(output_dir,"A"));
   rl= stod(read_sim_params(output_dir,"rl"));
   ru= stod(read_sim_params(output_dir,"ru"));
   exc_i = stoi(read_sim_params(output_dir, "exc_i"));
   collapse_and_bh = stoi(read_sim_params(output_dir, "collapse_and_bh"));

}
/*===========================================================================*/
Sim_params::~Sim_params(void)
{
}
