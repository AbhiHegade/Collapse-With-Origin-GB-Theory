#include <vector>
using std::vector;
#include <string>
using std::string;

#include "fd_stencils.hpp"



double Dx_ptpc_2nd(const double vp1, const double vm1, const double dx)
{

return (vp1 - vm1)/(2.*dx);

}

double Dx_ptp0_2nd(const double vp2, const double vp1, const double v0, const double dx)
{
  return (-vp2 + 4.*vp1 -3.*v0)/(2.*dx);
}

double Dx_2_ptpc_2nd(const double vp1 , const double v0, const double vm1, const double dx)
{
  return (vp1 - 2.*v0 + vm1)/(dx*dx);
}










// //=================================================================================
// inline double Dx_ptpc_2nd(const double vp1, const double vm1, const double dx);
// //=================================================================================
// inline double Dx_ptp0_2nd(const double vp2, const double vp1, const double v0, const double dx);
// //=================================================================================
// inline double Dx_2_ptpc_2nd(const double vp1, const double v0, const double vm1, const double dx);
/*===========================================================================*/
// inline double Dx_ptm1_4th(
//    const double vp1, const double v0, const double vm1,
//    const double vm2, const double vm3,
//    const double dx);
// /*===========================================================================*/
// inline double Dx_ptc_4th(
//    const double vp2, const double vp1, const double vm1, const double vm2,
//    const double dx);
// /*===========================================================================*/
// inline double Dx_ptp1_4th(
//    const double vp3, const double vp2, const double vp1,
//    const double v0, const double vm1,
//    const double dx);
// /*===========================================================================*/
// inline double Dx_ptp0_4th(
//    const double vp4, const double vp3, const double vp2,
//    const double vp1, const double v0,
//    const double dx);
// /*===========================================================================*/
// inline double make_Dx_zero(
//    const double vp4, const double vp3, const double vp2,
//    const double vp1);
// /*===========================================================================*/
// inline double Dx_2_ptpc_4th(
// 	const double vp3, const double vp2, const double vp1,
// 	const double vp0, const double vm1, const double vm2,
// 	const double vm3, const double dr
// );
//=================================================================================
