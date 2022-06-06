#include <vector>
using std::vector;
#include <string>
using std::string;

#include "fd_stencils.hpp"

/*===========================================================================*/
inline double Dx_ptm1_4th(
   const double vp1, const double v0, const double vm1,
   const double vm2, const double vm3,
   const double dx);
/*===========================================================================*/
inline double Dx_ptc_4th(
   const double vp2, const double vp1, const double vm1, const double vm2,
   const double dx);
/*===========================================================================*/
inline double Dx_ptp1_4th(
   const double vp3, const double vp2, const double vp1,
   const double v0, const double vm1,
   const double dx);
/*===========================================================================*/
inline double Dx_ptp0_4th(
   const double vp4, const double vp3, const double vp2,
   const double vp1, const double v0,
   const double dx);
/*===========================================================================*/
inline double make_Dx_zero(
   const double vp4, const double vp3, const double vp2,
   const double vp1);
/*===========================================================================*/
inline double Dx_2_ptpc_4th(
	const double vp3, const double vp2, const double vp1,
	const double vp0, const double vm1, const double vm2,
	const double vm3, const double dr
);
