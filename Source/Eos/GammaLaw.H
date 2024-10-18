#ifndef GAMMALAW_H
#define GAMMALAW_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "PhysicsConstants.H"

namespace pele::physics::eos {

struct GammaLaw
{
  using eos_type = GammaLaw;

  static std::string identifier() { return "GammaLaw"; }

  // ReallArrayLike can be anything with an [] operator that returns
  // amrex::Real&
  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
  RY2R(const RealArrayLike& RY, amrex::Real& R, int first = 0)
  {
    R = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      R += RY[first + i];
    }
  }

  // ReallArrayLike can be anything with an [] operator that returns
  // amrex::Real&
  template <typename RealArrayLike1, typename RealArrayLike2>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void RY2RRinvY(
    const RealArrayLike1& RY,
    amrex::Real& R,
    amrex::Real& Rinv,
    RealArrayLike2&& Y)
  {
    RY2R(RY, R);
    Rinv = 1.0 / R;
    for (int i = 0; i < NUM_SPECIES; i++) {
      Y[i] = RY[i] * Rinv;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RY2dRdY(
    const amrex::Real& /*R*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*dRdY[]*/)
  {
    amrex::Error("RY2dRdY not applicable for GammaLaw EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void molecular_weight(amrex::Real mw[NUM_SPECIES])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      mw[n] = Constants::AIRMW;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void inv_molecular_weight(amrex::Real imw[NUM_SPECIES])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      imw[n] = 1. / Constants::AIRMW;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void Y2WBAR(const amrex::Real* /*Y*/, amrex::Real& WBAR)
  {
    WBAR = Constants::AIRMW;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Ei(const amrex::Real T, amrex::Real Ei[NUM_SPECIES]) const
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
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real* /*Y[]*/,
    amrex::Real Ei[]) const
  {
    T2Ei(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(
    const amrex::Real /*T*/,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cv) const
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    Cv = Constants::RU / (wbar * (gamma - 1.0));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cv) const
  {
    TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cp) const
  {
    amrex::Real cv;
    TY2Cv(T, Y, cv);
    Cp = gamma * cv;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cp) const
  {
    TY2Cp(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Cpi(const amrex::Real /*T*/, amrex::Real Cpi[NUM_SPECIES]) const
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      Cpi[n] = gamma * Constants::RU / (Constants::AIRMW * (gamma - 1.0));
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(
    const amrex::Real R,
    const amrex::Real* /*Y*/,
    const amrex::Real E,
    const amrex::Real /*T*/,
    amrex::Real& P) const
  {
    P = (gamma - 1.0) * R * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPY2Cs(
    const amrex::Real R,
    const amrex::Real P,
    const amrex::Real* /*Y*/,
    amrex::Real& Cs) const
  {
    Cs = std::sqrt(gamma * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cs(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cs) const
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
  void EY2T(
    const amrex::Real E, const amrex::Real Y[NUM_SPECIES], amrex::Real& T) const
  {
    amrex::Real poverrho, wbar;
    poverrho = (gamma - 1.0) * E;
    Y2WBAR(Y, wbar);
    T = poverrho * wbar / Constants::RU;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(
    const amrex::Real /*R*/,
    const amrex::Real E,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& T) const
  {
    EY2T(E, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(const amrex::Real H, const amrex::Real* /*Y*/, amrex::Real& T) const
  {
    amrex::Real wbar = Constants::AIRMW;
    const amrex::Real Cv = Constants::RU / (wbar * (gamma - 1.0));
    T = H / (Cv * gamma);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(
    const amrex::Real /*R*/,
    const amrex::Real H,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& T) const
  {
    HY2T(H, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& P) const
  {
    amrex::Real Cv, E;
    TY2Cv(T, Y, Cv);
    E = Cv * T;
    P = (gamma - 1.0) * R * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RYP2T(
    const amrex::Real R,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real P,
    amrex::Real& T)
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    T = P * wbar / (R * Constants::RU);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2WDOT(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*WDOT[]*/)
  {
    amrex::Error("RTY2WDOT not physically possible for GammaLaw EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2JAC(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*JAC[]*/,
    const int /*HP*/)
  {
    amrex::Error("RTY2JAC not supported with GammaLaw EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2C(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*acti[]*/)
  {
    amrex::Error("RTY2C not supported with GammaLaw EoS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void Y2X(const amrex::Real Y[NUM_SPECIES], amrex::Real X[NUM_SPECIES])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      X[n] = Y[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void X2Y(const amrex::Real X[NUM_SPECIES], amrex::Real Y[NUM_SPECIES])
  {
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y[n] = X[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE(
    const amrex::Real P,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real T,
    amrex::Real& R,
    amrex::Real& E) const
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    R = P * wbar / (Constants::RU * T);
    E = P / (R * (gamma - 1.0));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void PYT2R(
    const amrex::Real P,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real T,
    amrex::Real& R)
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    R = P * wbar / (Constants::RU * T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2E(
    const amrex::Real R,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real P,
    amrex::Real& E) const
  {
    amrex::Real wbar;
    Y2WBAR(Y, wbar);
    E = P / (R * (gamma - 1.0));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2E(
    const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& E) const
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
  void RTY2E(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& E) const
  {
    TY2E(T, Y, E);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(const amrex::Real T, amrex::Real Hi[NUM_SPECIES]) const
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
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real* /*Y[]*/,
    amrex::Real Hi[NUM_SPECIES]) const
  {
    T2Hi(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void S(amrex::Real& s) { s = 1.0; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  TY2G(const amrex::Real /*T*/, const amrex::Real* /*Y*/, amrex::Real& G) const
  {
    G = gamma;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2G(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& G) const
  {
    TY2G(T, Y, G);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(
    const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H) const
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
  void RTY2H(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& H) const
  {
    TY2H(T, Y, H);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(
    const amrex::Real /*R*/,
    const amrex::Real /*P*/,
    const amrex::Real E,
    amrex::Real& dpdr_e) const
  {
    dpdr_e = (gamma - 1.0) * E;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RG2dpde(const amrex::Real R, const amrex::Real /*G*/, amrex::Real& dpde) const
  {
    dpde = (gamma - 1.0) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& dpde,
    amrex::Real& dpdr_e) const
  {
    amrex::Real G, P;
    TY2G(T, Y, G);
    RTY2P(R, T, Y, P);
    dpde = (G - 1.0) * R;
    dpdr_e = P / R;
  }

  AMREX_GPU_HOST_DEVICE
  GammaLaw(const EosParm<GammaLaw>* eparm)
  {
    // This should eventually error out if eparm is nullptr
    // But for now we require that eparm->gamma=gamma anyway,
    // and this supports some potential legacy use cases for now
    gamma = (eparm != nullptr) ? eparm->gamma : gamma;
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE explicit GammaLaw(Args... /*unused*/)
  {
  }

  amrex::Real gamma{Constants::gamma};
};

} // namespace pele::physics::eos
#endif
