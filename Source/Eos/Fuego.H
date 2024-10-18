#ifndef FUEGO_H
#define FUEGO_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include "mechanism.H"
#include "PhysicsConstants.H"

namespace pele::physics::eos {

struct Fuego
{
  using eos_type = Fuego;

  static std::string identifier() { return "Fuego"; }

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
    amrex::Error("RY2dRdY not applicable for Fuego EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void molecular_weight(amrex::Real mw[NUM_SPECIES]) { get_mw(mw); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void inv_molecular_weight(amrex::Real imw[NUM_SPECIES])
  {
    get_imw(imw);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2Cp(const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& Cp)
  {
    CKCPBS(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Cp(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cp)
  {
    TY2Cp(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2Cv(const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& Cv)
  {
    CKCVBS(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Cv(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cv)
  {
    TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void T2Cpi(const amrex::Real T, amrex::Real Cpi[NUM_SPECIES])
  {
    CKCPMS(T, Cpi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RPY2Cs(
    const amrex::Real R,
    const amrex::Real P,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cs)
  {
    amrex::Real tmp[NUM_SPECIES];
    amrex::Real wbar = 0.0;
    CKMMWY(Y, wbar);
    amrex::Real T = P * wbar / (R * Constants::RU);
    CKCVMS(T, tmp);
    amrex::Real Cv = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      Cv += Y[i] * tmp[i];
    }
    amrex::Real G = (wbar * Cv + Constants::RU) / (wbar * Cv);
    Cs = std::sqrt(G * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Cs(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& Cs)
  {
    amrex::Real tmp[NUM_SPECIES];
    amrex::Real P = 0.0;
    CKPY(R, T, Y, P);
    CKCVMS(T, tmp);
    amrex::Real Cv = 0.0;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      Cv += Y[i] * tmp[i];
    }
    amrex::Real wbar = 0.0;
    CKMMWY(Y, wbar);
    amrex::Real G = (wbar * Cv + Constants::RU) / (wbar * Cv);
    Cs = std::sqrt(G * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  EY2T(const amrex::Real E, const amrex::Real Y[NUM_SPECIES], amrex::Real& T)
  {
    // For Fuego this function is really just a wrapper for GET_T_GIVEN_EY
    // In SRK this will be different probably
    int lierr = 0;
    GET_T_GIVEN_EY(E, Y, T, lierr);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void REY2T(
    const amrex::Real /*R*/,
    const amrex::Real E,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& T)
  {
    EY2T(E, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  HY2T(const amrex::Real H, const amrex::Real Y[NUM_SPECIES], amrex::Real& T)
  {
    // For Fuego this function is really just a wrapper for GET_T_GIVEN_HY
    // In SRK this will be different probably
    int lierr = 0;
    GET_T_GIVEN_HY(H, Y, T, lierr);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RHY2T(
    const amrex::Real /*R*/,
    const amrex::Real H,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& T)
  {
    HY2T(H, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RYET2P(
    const amrex::Real R,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real E,
    amrex::Real& T,
    amrex::Real& P)
  {
    int lierr = 0;
    GET_T_GIVEN_EY(E, Y, T, lierr);
    CKPY(R, T, Y, P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2P(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& P)
  {
    CKPY(R, T, Y, P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RYP2T(
    const amrex::Real R,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real P,
    amrex::Real& T)
  {
    amrex::Real wbar = 0.0;
    CKMMWY(Y, wbar);
    T = P * wbar / (R * Constants::RU);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2WDOT(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real WDOT[NUM_SPECIES])
  {
    amrex::Real C[NUM_SPECIES];
    CKYTCR(R, T, Y, C);
    CKWC(T, C, WDOT);

    for (int n = 0; n < NUM_SPECIES; n++) {
      WDOT[n] *= mw(n);
    }
  }

#ifdef AMREX_USE_HIP
  // avoid register spillage on AMD, might be able to remove in future
  __attribute__((noinline))
#else
  AMREX_FORCE_INLINE
#endif
  AMREX_GPU_HOST_DEVICE static void
  RTY2JAC(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real JAC[],
    const int HP)
  {
    amrex::Real C[NUM_SPECIES];
    CKYTCR(R, T, Y, C);
    DWDOT(JAC, C, &T, &HP);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2C(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real acti[])
  {
    CKYTCR(R, T, Y, acti);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void T2Ei(const amrex::Real T, amrex::Real Ei[NUM_SPECIES])
  {
    CKUMS(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Ei(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real* /*Y[]*/,
    amrex::Real Ei[NUM_SPECIES])
  {
    T2Ei(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void Y2X(const amrex::Real Y[NUM_SPECIES], amrex::Real X[NUM_SPECIES])
  {
    CKYTX(Y, X);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void X2Y(const amrex::Real X[NUM_SPECIES], amrex::Real Y[NUM_SPECIES])
  {
    CKXTY(X, Y);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void PYT2RE(
    const amrex::Real P,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real T,
    amrex::Real& R,
    amrex::Real& E)
  {
    CKRHOY(P, T, Y, R);
    amrex::Real ei[NUM_SPECIES];
    CKUMS(T, ei);
    E = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      E += Y[n] * ei[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void PYT2R(
    const amrex::Real P,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real T,
    amrex::Real& R)
  {
    CKRHOY(P, T, Y, R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RYP2E(
    const amrex::Real R,
    const amrex::Real Y[NUM_SPECIES],
    const amrex::Real P,
    amrex::Real& E)
  {
    amrex::Real wbar = 0.0;
    CKMMWY(Y, wbar);
    amrex::Real T = P * wbar / (R * Constants::RU);
    amrex::Real ei[NUM_SPECIES];
    CKUMS(T, ei);
    E = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      E += Y[n] * ei[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2E(const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& E)
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
  static void RTY2E(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& E)
  {
    TY2E(T, Y, E);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void T2Hi(const amrex::Real T, amrex::Real Hi[NUM_SPECIES])
  {
    CKHMS(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Hi(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real* /*Y[]*/,
    amrex::Real Hi[NUM_SPECIES])
  {
    T2Hi(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void S(amrex::Real& s) { s = 1.0; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2G(const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& G)
  {
    amrex::Real wbar = 0.0, Cv, Cvx;
    TY2Cv(T, Y, Cv);
    CKMMWY(Y, wbar);
    Cvx = wbar * Cv;
    G = (Cvx + Constants::RU) / Cvx;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2G(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& G)
  {
    TY2G(T, Y, G);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2H(const amrex::Real T, const amrex::Real Y[NUM_SPECIES], amrex::Real& H)
  {
    amrex::Real hi[NUM_SPECIES];
    CKHMS(T, hi);
    H = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      H += Y[n] * hi[n];
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2H(
    const amrex::Real /*R*/,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& H)
  {
    TY2H(T, Y, H);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void Y2WBAR(const amrex::Real Y[NUM_SPECIES], amrex::Real& WBAR)
  {
    amrex::Real summ = 0.0;
    for (int i = 0; i < NUM_SPECIES; i++) {
      summ += Y[i] * imw(i);
    }
    WBAR = 1.0 / summ;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RPE2dpdr_e(
    const amrex::Real R,
    const amrex::Real P,
    const amrex::Real /* E */,
    amrex::Real& dpdr_e)
  {
    dpdr_e = P / R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  RG2dpde(const amrex::Real R, const amrex::Real G, amrex::Real& dpde)
  {
    dpde = (G - 1.0) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2dpde_dpdre(
    const amrex::Real R,
    const amrex::Real T,
    const amrex::Real Y[NUM_SPECIES],
    amrex::Real& dpde,
    amrex::Real& dpdr_e)
  {
    amrex::Real G, P;
    TY2G(T, Y, G);
    RTY2P(R, T, Y, P);
    dpde = (G - 1.0) * R;
    dpdr_e = P / R;
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE explicit Fuego(Args... /*unused*/)
  {
  }
};
} // namespace pele::physics::eos
#endif
