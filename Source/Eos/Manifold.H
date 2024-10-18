#ifndef MANIFOLD_H
#define MANIFOLD_H

#include "BlackBoxFunctionFactory.H"

namespace pele::physics::eos {

// Default - runtime selectable Table vs. Network for MANIFOLD EOS
// But if not using manifold, just keep 1D table for simplicity
#ifndef MANIFOLD_EOS_TYPE
#define MANIFOLD_EOS_TYPE 1
#endif
#ifndef MANIFOLD_DIM
#define MANIFOLD_DIM 1
#endif

// Compile time choice for Table vs. Network vs. Runtime Selectable
// Note right now Manifold EOS not supported for SYCL at all
// Runtime selection is  significant performnce hit
#if MANIFOLD_EOS_TYPE == 1
using ManifoldFunctionType = TabulatedFunction<MANIFOLD_DIM>;
#elif MANIFOLD_EOS_TYPE == 0
using ManifoldFunctionType = BlackBoxFunction;
#elif MANIFOLD_EOS_TYPE == 2
using ManifoldFunctionType = NeuralNetFunction;
#else
static_assert(false, "Invalid MANIFOLD_EOS_TYPE specified");
#endif

struct Manifold
{
  using eos_type = Manifold;

  static std::string identifier() { return "Manifold"; }

  // ReallArrayLike can be anything with an [] operator that returns
  // amrex::Real&
  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void
  RY2R(const RealArrayLike& RY, amrex::Real& R, int first = 0)
  {
    R = RY[first + NUM_SPECIES - 1];
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
  static void molecular_weight(amrex::Real* /*mw[]*/)
  {
    amrex::Error("molecular_weight for species does not have significance for "
                 "Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void inv_molecular_weight(amrex::Real* /*imw[]*/)
  {
    amrex::Error("inv_molecular_weight for species does not have significance "
                 "for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void Y2WBAR(const amrex::Real* /*Y*/, amrex::Real& /*WBAR*/)
  {
    amrex::Error("Y2WBAR is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Ei(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*Ei[]*/)
  {
    amrex::Error("RTY2Ei does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2E(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real /*E*/)
  {
    amrex::Error("RTY2E does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void T2Hi(const amrex::Real /*T*/, amrex::Real Hi[])
  {
    // TODO: FIXME (need to remove enthalpy calculations from
    // PeleLM)
    for (int i = 0; i < NUM_SPECIES; i++) {
      Hi[i] = 0.0;
    }
    // amrex::Error("T2Hi does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Hi(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*Hi[]*/)
  {
    amrex::Error("RTY2Hi does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void TY2Cv(
    const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& /*Cv*/)
  {
    amrex::Error("TY2Cv is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Cv(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real& /*Cv*/)
  {
    amrex::Error("RTY2Cv is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2Cp(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& Cp)
  {
    // TODO: FIXME (need to remove enthalpy calculations from
    // PeleLM)
    Cp = 100.0;
    // amrex::Error("TY2Cp is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2Cp(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real& Cp)
  {
    // TODO: FIXME (need to remove enthalpy calculations from
    // PeleLM)
    Cp = 100.0;
    // amrex::Error("TY2Cp is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(
    const amrex::Real /*R*/,
    const amrex::Real /*E*/,
    const amrex::Real Y[],
    amrex::Real& T)
  {
    if (eosparm->compute_temperature) {
      manfunc.get_func()->get_value(eosparm->idx_T, Y, T);
    } else {
      T = 1.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(const amrex::Real /*H*/, const amrex::Real Y[], amrex::Real& T)
  {
    if (eosparm->compute_temperature) {
      manfunc.get_func()->get_value(eosparm->idx_T, Y, T);
    } else {
      T = 1.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(
    const amrex::Real /*R*/,
    const amrex::Real /*H*/,
    const amrex::Real Y[],
    amrex::Real& T)
  {
    if (eosparm->compute_temperature) {
      manfunc.get_func()->get_value(eosparm->idx_T, Y, T);
    } else {
      T = 1.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  TY2H(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& H)
  {
    H = 0.0;
    // TODO: FIXME Eliminate this call from PeleLM?
    // amrex::Error("TY2H is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void X2Y(const amrex::Real X[], amrex::Real Y[])
  {
    // TODO: FIXME
    // just return Y as a hack for PMF, this should really not be needed
    // anywhere
    for (int i = 0; i < NUM_SPECIES; i++) {
      Y[i] = X[i];
    }
    // amrex::Error("X2Y does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void Y2X(const amrex::Real* /*Y[]*/, amrex::Real* /*X[]*/)
  {
    amrex::Error("Y2X does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2WDOT(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real Y[],
    amrex::Real WDOT[])
  {
    manfunc.get_func()->get_values(NUM_SPECIES - 1, eosparm->idx_Wdot, Y, WDOT);

    // Density source is 0
    WDOT[NUM_SPECIES - 1] = 0.0;
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
    amrex::Error("RTY2JAC is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RTY2C(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*JAC[]*/)
  {
    amrex::Error("RTY2C does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(
    const amrex::Real R,
    const amrex::Real /*T*/,
    const amrex::Real Y[],
    amrex::Real& P)
  {
    // Get value of Rho consistent with other state variables
    amrex::Real rho_from_table;
    manfunc.get_func()->get_value(eosparm->idx_density, Y, rho_from_table);
    switch (eosparm->dens_lookup) {
    case density_lookup_type::linear:
      break;
    case density_lookup_type::log:
      rho_from_table = exp(rho_from_table);
      break;
    case density_lookup_type::inverse:
      rho_from_table = 1.0 / rho_from_table;
      break;
    }

    // Return pressure is nominal pressure scaled due to density mismatch
    P = eosparm->Pnom_cgs * R / rho_from_table;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(
    const amrex::Real /*P*/,
    const amrex::Real Y[],
    const amrex::Real /*T*/,
    amrex::Real& R)
  {
    manfunc.get_func()->get_value(eosparm->idx_density, Y, R);
    switch (eosparm->dens_lookup) {
    case density_lookup_type::linear:
      break;
    case density_lookup_type::log:
      R = exp(R);
      break;
    case density_lookup_type::inverse:
      R = 1.0 / R;
      break;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RY2dRdY(const amrex::Real R, const amrex::Real Y[], amrex::Real dRdY[])
  {
    manfunc.get_func()->get_derivs(eosparm->idx_density, Y, dRdY);
    switch (eosparm->dens_lookup) {
    case density_lookup_type::linear:
      break;
    case density_lookup_type::log:
      for (int i = 0; i < NUM_SPECIES - 1; i++) {
        dRdY[i] *= R;
      }
      break;
    case density_lookup_type::inverse:
      for (int i = 0; i < NUM_SPECIES - 1; i++) {
        dRdY[i] *= -R * R;
      }
      break;
    }
  }

  AMREX_GPU_HOST_DEVICE
  Manifold(const EosParm<Manifold>* eparm)
    : eosparm(eparm), manfunc(eparm->manf_data)
  {
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE Manifold()
  {
    amrex::Error("Manifold EOS must be initialized with table data or neural "
                 "network via EosParm");
  }

  AMREX_GPU_HOST_DEVICE
  ~Manifold() = default;

private:
  const EosParm<Manifold>* eosparm;
  BlackBoxFunctionFactory<ManifoldFunctionType> manfunc;
};

} // namespace pele::physics::eos
#endif
