#ifndef EOS_H_
#define EOS_H_

#include <AMReX.H>
#include <AMReX_REAL.H> 
#include <AMReX_Utility.H>
#include "mechanism.h"
#define naux 0
 
using namespace amrex;
 
class EOS
{
public:
  AMREX_GPU_HOST_DEVICE EOS();
  AMREX_GPU_HOST_DEVICE ~EOS();

  AMREX_GPU_HOST_DEVICE void eos_EY2T(const amrex::Real *Y, amrex::Real &E, amrex::Real &T);
  AMREX_GPU_HOST_DEVICE void eos_T2EI(      amrex::Real  T, amrex::Real *ei);
  AMREX_GPU_HOST_DEVICE void eos_TY2Cv(     amrex::Real  T, const amrex::Real *Y, amrex::Real &Cv);

  AMREX_GPU_HOST_DEVICE void eos_HY2T(const amrex::Real *Y,       amrex::Real  H, amrex::Real &T);
  AMREX_GPU_HOST_DEVICE void eos_T2HI(      amrex::Real  T,       amrex::Real *hi);
  AMREX_GPU_HOST_DEVICE void eos_TY2Cp(     amrex::Real  T, const amrex::Real *Y, amrex::Real &Cp);

  AMREX_GPU_HOST_DEVICE void eos_RTY2C(  amrex::Real rho, amrex::Real T, const amrex::Real *Y, amrex::Real *acti);
  AMREX_GPU_HOST_DEVICE void eos_RTY2W(  amrex::Real rho, amrex::Real T, const amrex::Real *Y, amrex::Real *wdot);
  AMREX_GPU_HOST_DEVICE void eos_RTY2JAC(amrex::Real rho, amrex::Real T, const amrex::Real *Y, amrex::Real *Jac,  int HP);
};

#endif
