#include "Fuego_EOS.H"
#include "mechanism.h"
#include <chemistry_file.H>
#include <cmath>



AMREX_GPU_HOST_DEVICE EOS::EOS()
{}

AMREX_GPU_HOST_DEVICE EOS::~EOS()
{}


AMREX_GPU_HOST_DEVICE
void eos_EY2T(const Real *Y, Real E, Real &T)
{
  int lierr=0;
  GET_T_GIVEN_EY(&E, (Real*)Y, &T, &lierr);
}

AMREX_GPU_HOST_DEVICE
void eos_T2EI(Real T, Real *ei)
{
  CKUMS(&T,  ei);
}

AMREX_GPU_HOST_DEVICE
void eos_TY2Cv(Real T, const Real *Y, Real &Cv)
{
  Real cvi[NUM_SPECIES];

  CKCVMS(&T,  cvi);

  Cv = 0.0;
  for(int i = 0; i < NUM_SPECIES; ++i){
    Cv += Y[i] * cvi[i];
  }
}


AMREX_GPU_HOST_DEVICE
void eos_HY2T(const Real *Y, Real H, Real &T)
{
  int lierr=0;
  GET_T_GIVEN_HY(&H, (Real *)Y, &T, &lierr);
}

AMREX_GPU_HOST_DEVICE
void eos_T2HI(Real T, Real *hi)
{
  CKHMS(&T,  hi);
}

AMREX_GPU_HOST_DEVICE
void eos_TY2Cp(Real T, const Real *Y, Real &Cp)
{
  Real cpi[NUM_SPECIES];

  CKCPMS(&T,  cpi);

  Cp = 0.0;
  for(int i = 0; i < NUM_SPECIES; ++i){
    Cp += Y[i] * cpi[i];
  }
}

AMREX_GPU_HOST_DEVICE
void eos_RTY2C(Real rho, Real T, const Real *Y, Real *acti)
{
  CKYTCR(&rho, &T, (Real *)Y, acti);
}

AMREX_GPU_HOST_DEVICE
void eos_RTY2W(Real rho, Real T, const Real *Y, Real *wdot)
{
  Real C[NUM_SPECIES];

  CKYTCR(&rho, &T, (Real *)Y, C);
  CKWC(&T, C, wdot);
}

/* Should not be here but I am not sure wether the CKYTCR needs to be wraped. */
AMREX_GPU_HOST_DEVICE
void eos_RTY2JAC(Real rho, Real T, const Real *Y, Real *Jac,  int HP)
{
  Real C[NUM_SPECIES];

  CKYTCR(&rho, &T, (Real *)Y, C);
  DWDOT(Jac,C,&T,&HP);
}
