#ifndef _EOS_H_
#define _EOS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include "chemistry_file.H"

namespace EOS {

extern AMREX_GPU_DEVICE_MANAGED amrex::Real gamma;
constexpr amrex::Real RU = 8.31446261815324e7;
constexpr amrex::Real RUC  = 1.98721558317399615845;
constexpr amrex::Real PATM = 1.01325e+06;
constexpr amrex::Real AIRMW = 28.97;

void init();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
molecular_weight(amrex::Real mw[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    mw[n] = AIRMW;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
inv_molecular_weight(amrex::Real imw[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    imw[n] = 1./AIRMW;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Y2WBAR(amrex::Real Y[], amrex::Real& WBAR)
{
  WBAR = AIRMW;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Ei(amrex::Real T, amrex::Real Ei[])
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  for (int n = 0; n < NUM_SPECIES; n++)
    Ei[n] = Cv * T;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2Cv(amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
{
  amrex::Real wbar;
  EOS::Y2WBAR(Y, wbar);
  Cv = RU / (wbar * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
{
  amrex::Real cv;
  EOS::TY2Cv(T, Y, cv);
  Cp = gamma * cv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Cpi(amrex::Real T, amrex::Real Cpi[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    Cpi[n] = gamma * RU / (AIRMW * (gamma - 1.0));
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
  P = (gamma - 1.0) * R * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real Y[], amrex::Real& Cs)
{
  Cs = std::sqrt(gamma * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs)
{
  amrex::Real wbar, P;
  EOS::Y2WBAR(Y, wbar);
  amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  amrex::Real E = Cv * T;
  EOS::RYET2P(R, Y, E, T, P);
  Cs = std::sqrt(gamma * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T)
{
  amrex::Real poverrho, wbar;
  poverrho = (gamma - 1.0) * E;
  EOS::Y2WBAR(Y, wbar);
  T = poverrho * wbar / RU;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HY2T(amrex::Real H, amrex::Real Y[], amrex::Real& T)
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  T = H / (Cv * gamma);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TRY2P(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& P)
{
  amrex::Real Cv, E;
  EOS::TY2Cv(T, Y, Cv);
  E = Cv * T;
  P = (gamma - 1.0) * R * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Y2X(amrex::Real Y[], amrex::Real X[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    X[n] = Y[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
X2Y(amrex::Real X[], amrex::Real Y[])
{
  for (int n = 0; n < NUM_SPECIES; n++)
    Y[n] = X[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
PYTR2E(
  amrex::Real P, amrex::Real Y[], amrex::Real T, amrex::Real& R, amrex::Real& E)
{
  amrex::Real wbar;
  EOS::Y2WBAR(Y, wbar);
  R = P * wbar / (RU * T);
  E = P / (R * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RYP2E(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& E)
{
  amrex::Real wbar;
  EOS::Y2WBAR(Y, wbar);
  E = P / (R * (gamma - 1.0));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Hi(amrex::Real T, amrex::Real Hi[])
{
  amrex::Real wbar = AIRMW;
  const amrex::Real Cv = RU / (wbar * (gamma - 1.0));
  for (int n = 0; n < NUM_SPECIES; n++)
    Hi[n] = Cv * T * gamma;
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
  G = gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2G(amrex::Real T, amrex::Real Y[NUM_SPECIES], amrex::Real& G)
{
  G = gamma;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPE2dpdr_e(amrex::Real R, amrex::Real P, amrex::Real E, amrex::Real& dpdr_e)
{
  dpdr_e = (gamma - 1.0) * E;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RG2dpde(amrex::Real R, amrex::Real G, amrex::Real& dpde)
{
  dpde = (gamma - 1.0) * R;
}

}; // namespace EOS

#endif
