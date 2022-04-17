#ifndef _GAMMALAW_H_
#define _GAMMALAW_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "PhysicsConstants.H"

namespace pele {
namespace physics {
namespace eos {

struct GammaLaw
{
  using eos_type = GammaLaw;

  static std::string identifier() { return "GammaLaw"; }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void molecular_weight(amrex::Real mw[])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      mw[n] = Constants::AIRMW;
    }
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void inv_molecular_weight(amrex::Real imw[])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      imw[n] = 1. / Constants::AIRMW;
    }
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void Y2WBAR(amrex::Real* /*Y*/, amrex::Real& WBAR)
  {
    WBAR = Constants::AIRMW;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[]) const
  {
    amrex::Real wbar = Constants::AIRMW;
    const amrex::Real Cv = Constants::RU / (wbar * (gamma - 1.0));
    for (int n = 0; n < NUM_SPECIES; n++) {
      Ei[n] = Cv * T;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(
    amrex::Real /*R*/,
    amrex::Real T,
    amrex::Real* /*Y[]*/,
    amrex::Real Ei[]) const
  {
    T2Ei(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real /*T*/, amrex::Real Y[], amrex::Real& Cv) const
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    Cv = Constants::RU / (wbar * (gamma - 1.0));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(
    amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& Cv) const
  {
    TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp) const
  {
    amrex::Real cv;
    TY2Cv(T, Y, cv);
    Cp = gamma * cv;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(
    amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& Cp) const
  {
    TY2Cp(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real /*T*/, amrex::Real Cpi[]) const
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      Cpi[n] = gamma * Constants::RU / (Constants::AIRMW * (gamma - 1.0));
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(
    amrex::Real R,
    amrex::Real* /*Y*/,
    amrex::Real& E,
    amrex::Real& /*T*/,
    amrex::Real& P) const
  {
    P = (gamma - 1.0) * R * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPY2Cs(
    amrex::Real R, amrex::Real P, amrex::Real* /*Y*/, amrex::Real& Cs) const
  {
    Cs = std::sqrt(gamma * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs) const
  {
    amrex::Real wbar, P;
    Y2WBAR(Y, wbar);
    amrex::Real Cv = Constants::RU / (wbar * (gamma - 1.0));
    amrex::Real E = Cv * T;
    RYET2P(R, Y, E, T, P);
    Cs = std::sqrt(gamma * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T) const
  {
    amrex::Real poverrho, wbar;
    poverrho = (gamma - 1.0) * E;
    Y2WBAR(Y, wbar);
    T = poverrho * wbar / Constants::RU;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  REY2T(amrex::Real /*R*/, amrex::Real E, amrex::Real Y[], amrex::Real& T) const
  {
    EY2T(E, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real* /*Y*/, amrex::Real& T) const
  {
    amrex::Real wbar = Constants::AIRMW;
    const amrex::Real Cv = Constants::RU / (wbar * (gamma - 1.0));
    T = H / (Cv * gamma);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RHY2T(amrex::Real /*R*/, amrex::Real H, amrex::Real Y[], amrex::Real& T) const
  {
    HY2T(H, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2P(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& P) const
  {
    amrex::Real Cv, E;
    TY2Cv(T, Y, Cv);
    E = Cv * T;
    P = (gamma - 1.0) * R * E;
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void
  RYP2T(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& T)
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    T = P * wbar / (R * Constants::RU);
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void RTY2WDOT(
    amrex::Real /*R*/,
    amrex::Real /*T*/,
    amrex::Real* /*Y[]*/,
    amrex::Real* /*WDOT[]*/)
  {
    amrex::Error("RTY2WDOT not physically possible for GammaLaw EoS");
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void RTY2JAC(
    amrex::Real /*R*/,
    amrex::Real /*T*/,
    amrex::Real* /*Y[]*/,
    amrex::Real* /*JAC[]*/,
    const int /*HP*/)
  {
    amrex::Error("RTY2JAC not supported with GammaLaw EoS");
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void RTY2C(
    amrex::Real /*R*/,
    amrex::Real /*T*/,
    amrex::Real* /*Y[]*/,
    amrex::Real* /*acti[]*/)
  {
    amrex::Error("RTY2C not supported with GammaLaw EoS");
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void Y2X(const amrex::Real Y[], amrex::Real X[])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      X[n] = Y[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void X2Y(const amrex::Real X[], amrex::Real Y[])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y[n] = X[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE(
    amrex::Real P,
    amrex::Real Y[],
    amrex::Real T,
    amrex::Real& R,
    amrex::Real& E) const
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    R = P * wbar / (Constants::RU * T);
    E = P / (R * (gamma - 1.0));
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void
  PYT2R(amrex::Real P, amrex::Real Y[], amrex::Real T, amrex::Real& R)
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    R = P * wbar / (Constants::RU * T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RYP2E(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& E) const
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    E = P / (R * (gamma - 1.0));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2E(amrex::Real T, const amrex::Real Y[], amrex::Real& E) const
  {
    amrex::Real ei[NUM_SPECIES];
    T2Ei(T, ei);
    E = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      E += Y[n] * ei[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2E(amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& E) const
  {
    TY2E(T, Y, E);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[]) const
  {
    amrex::Real wbar = Constants::AIRMW;
    const amrex::Real Cv = Constants::RU / (wbar * (gamma - 1.0));
    for (int n = 0; n < NUM_SPECIES; n++) {
      Hi[n] = Cv * T * gamma;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(
    amrex::Real /*R*/,
    amrex::Real T,
    amrex::Real* /*Y[]*/,
    amrex::Real Hi[]) const
  {
    T2Hi(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  static AMREX_FORCE_INLINE void S(amrex::Real& s) { s = 1.0; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2G(amrex::Real /*T*/, amrex::Real* /*Y*/, amrex::Real& G) const
  {
    G = gamma;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2G(amrex::Real /*R*/, amrex::Real T, amrex::Real* Y, amrex::Real& G) const
  {
    TY2G(T, Y, G);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  TY2H(amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H) const
  {
    amrex::Real Hi[NUM_SPECIES];
    amrex::Real wbar = Constants::AIRMW;
    const amrex::Real Cv = Constants::RU / (wbar * (gamma - 1.0));
    H = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      Hi[n] = Cv * T * gamma;
      H = H + Y[n] * Hi[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(
    amrex::Real /*R*/,
    amrex::Real /*P*/,
    amrex::Real E,
    amrex::Real& dpdr_e) const
  {
    dpdr_e = (gamma - 1.0) * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real /*G*/, amrex::Real& dpde) const
  {
    dpde = (gamma - 1.0) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(
    amrex::Real R,
    amrex::Real T,
    amrex::Real Y[],
    amrex::Real& dpde,
    amrex::Real& dpdr_e) const
  {
    amrex::Real G, P;
    TY2G(T, Y, G);
    RTY2P(R, T, Y, P);
    dpde = (G - 1.0) * R;
    dpdr_e = P / R;
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE GammaLaw(Args... /*unused*/)
  {
  }
  amrex::Real gamma{Constants::gamma};
};

} // namespace eos
} // namespace physics
} // namespace pele
#endif
