#ifndef _EOS_H_
#define _EOS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include "mechanism.h"
#include "chemistry_file.H"

//#include "IndexDefines.H"
//#include "Constants.H"

namespace EOS {

//extern AMREX_GPU_DEVICE_MANAGED int upassMap[NPASSIVE];
//extern AMREX_GPU_DEVICE_MANAGED int qpassMap[NPASSIVE];

void init();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
GET_MW(amrex::Real mw[])
{
  get_mw(mw);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
GET_IMW(amrex::Real imw[])
{
  get_imw(imw);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
{
  amrex::Real tmp[NUM_SPECIES];
  CKCPMS(&T, tmp);
  Cp = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i)
    Cp += Y[i] * tmp[i];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2Cv(amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
{
  amrex::Real tmp[NUM_SPECIES];
  CKCVMS(&T, tmp);
  Cv = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i)
    Cv += Y[i] * tmp[i];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Cpi(amrex::Real T, amrex::Real Cpi[])
{
  CKCPMS(&T, Cpi);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real Y[], amrex::Real& Cs)
{
  amrex::Real tmp[NUM_SPECIES];
  amrex::Real wbar;
  amrex::Real RU, RUC, PA;
  CKMMWY(&Y[0], &wbar);
  CKRP(&RU,&RUC,&PA);
  amrex::Real T = P * wbar / (R * RU);
  CKCVMS(&T, tmp);
  amrex::Real Cv = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++)
    Cv += Y[i] * tmp[i];
  amrex::Real G = (wbar * Cv + RU) / (wbar * Cv);
  Cs = std::sqrt(G * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs)
{
  amrex::Real tmp[NUM_SPECIES];
  amrex::Real RU, RUC, PA;
  amrex::Real P;
  CKPY(&R, &T, Y, &P);
  CKRP(&RU,&RUC,&PA);
  CKCVMS(&T, tmp);
  amrex::Real Cv = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i)
    Cv += Y[i] * tmp[i];
  amrex::Real wbar;
  CKMMWY(&Y[0], &wbar);
  amrex::Real G = (wbar * Cv + RU) / (wbar * Cv);
  Cs = std::sqrt(G * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T)
{
  // For Fuego this function is really just a wrapper for GET_T_GIVEN_EY
  // In SRK this will be different probably
  int lierr = 0;
  GET_T_GIVEN_EY(&E, Y, &T, &lierr);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HY2T(amrex::Real H, amrex::Real Y[], amrex::Real& T)
{
  // For Fuego this function is really just a wrapper for GET_T_GIVEN_HY
  // In SRK this will be different probably
  int lierr = 0;
  GET_T_GIVEN_HY(&H, Y, &T, &lierr);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RYET2P(
  amrex::Real R,
  amrex::Real Y[],
  amrex::Real& E,
  amrex::Real& T,
  amrex::Real& P)
{
  int lierr = 0;
  GET_T_GIVEN_EY(&E, Y, &T, &lierr);
  CKPY(&R, &T, Y, &P);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TRY2P(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& P)
{
  CKPY(&R, &T, Y, &P);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2WDOT(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real WDOT[])
{
  //CKWYR(&R, &T, Y, WDOT);
  amrex::Real C[NUM_SPECIES];
  CKYTCR(&R, &T, Y, C);
  CKWC(&T, C, WDOT);

  amrex::Real mw[NUM_SPECIES];
  get_mw(mw);
  for (int n = 0; n < NUM_SPECIES; n++)
    WDOT[n] *= mw[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2JAC(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real JAC[], int HP)
{
  amrex::Real C[NUM_SPECIES];
  CKYTCR(&R, &T, Y, C);
  DWDOT(JAC, C, &T,&HP);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2C(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real acti[])
{
  CKYTCR(&R, &T,  Y, acti); 
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
YT2E(amrex::Real Y[], amrex::Real T, amrex::Real E[])
{
  CKUMS(&T, E);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
YT2Cv(amrex::Real Y[], amrex::Real T, amrex::Real& Cv)
{
  Cv = 0.0;
  amrex::Real tmp[NUM_SPECIES];
  CKCVMS(&T, tmp);
  for (int i = 0; i < NUM_SPECIES; i++)
    Cv += Y[i] * tmp[i];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Y2X(amrex::Real Y[], amrex::Real X[])
{
  CKYTX(Y, X);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
X2Y(amrex::Real X[], amrex::Real Y[])
{
  CKXTY(X, Y);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
PYTR2E(
  amrex::Real P, amrex::Real Y[], amrex::Real T, amrex::Real& R, amrex::Real& E)
{
  CKRHOY(&P, &T, Y, &R);
  amrex::Real ei[NUM_SPECIES];
  CKUMS(&T, ei);
  E = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++)
    E += Y[n] * ei[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RYP2E(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& E)
{
  amrex::Real wbar;
  amrex::Real RU, RUC, PA;
  CKRP(&RU,&RUC,&PA);
  CKMMWY(&Y[0], &wbar);
  amrex::Real T = P * wbar / (R * RU);
  amrex::Real ei[NUM_SPECIES];
  CKUMS(&T, ei);
  E = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++)
    E += Y[n] * ei[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
YT2H(amrex::Real Y[], amrex::Real T, amrex::Real H[])
{
  CKHMS(&T, H);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
S(amrex::Real& s)
{
  s = 1.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
EY2G(amrex::Real E, amrex::Real Y[NUM_SPECIES], amrex::Real& G)
{
  amrex::Real wbar, Cv, T, Cvx;
  amrex::Real RU, RUC, PA;
  CKRP(&RU,&RUC,&PA);
  EOS::EY2T(E, Y, T);
  EOS::YT2Cv(Y, T, Cv);
  CKMMWY(&Y[0], &wbar);
  Cvx = wbar * Cv;
  G = (Cvx + RU) / Cvx;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Y2WBAR(amrex::Real Y[], amrex::Real& WBAR)
{
  amrex::Real tmp[NUM_SPECIES];
  get_imw(tmp);
  amrex::Real summ = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++)
    summ += Y[i] * tmp[i];
  WBAR = 1.0 / summ;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPE2dpdr_e(amrex::Real R, amrex::Real P, amrex::Real E, amrex::Real& dpdr_e)
{
  dpdr_e = P / R;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RG2dpde(amrex::Real R, amrex::Real G, amrex::Real& dpde)
{
  dpde = (G - 1.0) * R;
}

// Hydro -> Advected -> Species -> Aux
// If num_adv == 0 -> QFA = QFS and UFA = UFS, see IndexDefines.H
// For explicit definitions.
//AMREX_GPU_HOST_DEVICE
//AMREX_FORCE_INLINE
//int
//upass_map(const int i)
//{
//  return upassMap[i];
//}
//
//AMREX_GPU_HOST_DEVICE
//AMREX_FORCE_INLINE
//int
//qpass_map(const int i)
//{
//  /*V and W are passive*/
//  return qpassMap[i];
//}

}; // namespace EOS

#endif
