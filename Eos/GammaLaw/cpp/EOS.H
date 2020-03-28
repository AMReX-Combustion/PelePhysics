#ifndef _EOS_H_
#define _EOS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include "chemistry_file.H"

#include "IndexDefines.H"
#include "Constants.H"

namespace EOS {

extern AMREX_GPU_DEVICE_MANAGED int upassMap[NPASSIVE];
extern AMREX_GPU_DEVICE_MANAGED int qpassMap[NPASSIVE];
extern AMREX_GPU_DEVICE_MANAGED amrex::Real gamma;

void init();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_wb
Y2WBAR(amrex::Real Y[], amrex::Real& WBAR)
{
  const amrex::Real molec_wt[NUM_SPECIES] = {28.97};
  amrex::Real summ = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i)
    summ += Y[i] * (1.0 / molec_wt[i]);
  WBAR = 1.0 / summ;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_cv
YT2Cv(amrex::Real Y[], amrex::Real T, amrex::Real& Cv)
{
  amrex::Real wbar;
  Y2WBAR(Y, wbar);
  Cv = RU / (wbar * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_cp
TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
{
  amrex::Real cv;
  YT2Cv(Y, T, cv);
  Cp = gamma * cv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_cpi
T2Cpi(amrex::Real T, amrex::Real Cpi[])
{
  Cpi[0] = gamma*RU/ (28.97 * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_p
RYET2P(
  amrex::Real R,
  amrex::Real Y[],
  amrex::Real& E,
  amrex::Real& T,
  amrex::Real& P)
{
  P = (gamma - 1.0) * R * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_cs_no_T
RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real Y[], amrex::Real& Cs)
{
  Cs = std::sqrt(gamma * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_cs
RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs)
{
  amrex::Real wbar, P;
  Y2WBAR(Y, wbar);
  amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  amrex::Real E = Cv * T;
  RYET2P(R, Y, E, T, P);
  Cs = std::sqrt(gamma * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //cmpT
EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T)
{
  amrex::Real poverrho, wbar;
  poverrho = (gamma - 1.0) * E;
  Y2WBAR(Y, wbar);
  T = poverrho * wbar / RU;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //TrhoY2p
TRY2P(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& P)
{
  amrex::Real Cv, E;
  YT2Cv(Y, T, Cv);
  E = Cv * T;
  P = (gamma - 1.0) * R * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //ytx
Y2X(amrex::Real Y[], amrex::Real X[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    X[n] = Y[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //tp
PYTR2E(
  amrex::Real P,
  amrex::Real Y[],
  amrex::Real T,
  amrex::Real& R,
  amrex::Real& E)
{
  amrex::Real wbar;
  Y2WBAR(Y, wbar);
  R = P * wbar / (RU * T);
  E = P / (R * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //rhopY2e
RYP2E(
  amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& E)
{
  amrex::Real wbar;
  Y2WBAR(Y, wbar);
  E = P / (R * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_hi
YT2H(amrex::Real Y[], amrex::Real T, amrex::Real H[])
{
  amrex::Real wbar;
  Y2WBAR(Y, wbar);
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  for (int n = 0; n < NUM_SPECIES; n++)
    H[n] = Cv * T * gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_s
S(amrex::Real& s)
{
  s = 1.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void //get_gamma
EY2G(amrex::Real E, amrex::Real Y[NUM_SPECIES], amrex::Real& G)
{
  G = gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPE2dpdr_e(
  amrex::Real R,
  amrex::Real P,
  amrex::Real E,
  amrex::Real& dpdr_e)
{
  dpdr_e = (gamma - 1.0) * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RG2dpde(
  amrex::Real R,
  amrex::Real G,
  amrex::Real& dpde)
{
  dpde = (gamma - 1.0) * R;
}

// Hydro -> Advected -> Species -> Aux
// If num_adv == 0 -> QFA = QFS and UFA = UFS, see IndexDefines.H
// For explicit definitions.
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
upass_map(const int i)
{
  return upassMap[i];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
qpass_map(const int i)
{
  return qpassMap[i];
}

}; // namespace EOS

#endif
