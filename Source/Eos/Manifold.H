#ifndef _MANIFOLD_H_
#define _MANIFOLD_H_

#include <memory>
#include "Table.H"

#ifdef USE_LIBTORCH
#include "NeuralNet.H"
#else
#include "NeuralNetHomerolled.H"
#endif

namespace pele {
namespace physics {
namespace eos {

struct Manifold
{
  using eos_type = Manifold;

  static std::string identifier() { return "Manifold"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RY2R(const amrex::Real RY[], amrex::Real& R)
  {
    R = RY[NUM_SPECIES - 1];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RY2RRinvY(
    const amrex::Real RY[], amrex::Real& R, amrex::Real& Rinv, amrex::Real Y[])
  {
    RY2R(RY, R);
    Rinv = 1.0 / R;
    for (int is = 0; is < NUM_SPECIES; is++) {
      Y[is] = RY[is] * Rinv;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void molecular_weight(amrex::Real* /*mw[]*/)
  {
    amrex::Error("molecular_weight for species does not have significance for "
                 "Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void inv_molecular_weight(amrex::Real* /*imw[]*/)
  {
    amrex::Error("inv_molecular_weight for species does not have significance "
                 "for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(const amrex::Real* /*Y*/, amrex::Real& /*WBAR*/)
  {
    amrex::Error("Y2WBAR is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*Ei[]*/)
  {
    amrex::Error("RTY2Ei does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2E(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*E*/)
  {
    amrex::Error("RTY2E does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(const amrex::Real /*T*/, amrex::Real Hi[])
  {
    // TODO: FIXME (need to remove enthalpy calculations from
    // PeleLM)
    for (int is = 0; is < NUM_SPECIES; is++) {
      Hi[is] = 0.0;
    }
    // amrex::Error("T2Hi does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real* /*Hi[]*/)
  {
    amrex::Error("RTY2Hi does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(
    const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& /*Cv*/)
  {
    amrex::Error("TY2Cv is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(
    const amrex::Real /*R*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*Y[]*/,
    amrex::Real& /*Cv*/)
  {
    amrex::Error("RTY2Cv is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  TY2Cp(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& Cp)
  {
    // TODO: FIXME (need to remove enthalpy calculations from
    // PeleLM)
    Cp = 100.0;
    // amrex::Error("TY2Cp is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(
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
      manfunc->get_value(eosparm->idx_T, Y, T);
    } else {
      T = 1.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(const amrex::Real /*H*/, const amrex::Real Y[], amrex::Real& T)
  {
    if (eosparm->compute_temperature) {
      manfunc->get_value(eosparm->idx_T, Y, T);
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
      manfunc->get_value(eosparm->idx_T, Y, T);
    } else {
      T = 1.0;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& H)
  {
    H = 0.0;
    // TODO: FIXME Eliminate this call from PeleLM?
    // amrex::Error("TY2H is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void X2Y(const amrex::Real X[], amrex::Real Y[])
  {
    // TODO: FIXME
    // just return Y as a hack for PMF, this should really not be needed
    // anywhere
    for (int is = 0; is < NUM_SPECIES; is++) {
      Y[is] = X[is];
    }
    // amrex::Error("X2Y does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2X(const amrex::Real* /*Y[]*/, amrex::Real* /*X[]*/)
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
    manfunc->get_values(NUM_SPECIES - 1, eosparm->idx_Wdot, Y, WDOT);

    // Density source is 0
    WDOT[NUM_SPECIES - 1] = 0.0;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(
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
  void RTY2C(
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
    manfunc->get_value(eosparm->idx_density, Y, rho_from_table);
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
    manfunc->get_value(eosparm->idx_density, Y, R);
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
    manfunc->get_derivs(eosparm->idx_density, Y, dRdY);
    switch (eosparm->dens_lookup) {
    case density_lookup_type::linear:
      break;
    case density_lookup_type::log:
      for (int is = 0; is < NUM_SPECIES - 1; is++) {
        dRdY[is] *= R;
      }
      break;
    case density_lookup_type::inverse:
      for (int is = 0; is < NUM_SPECIES - 1; is++) {
        dRdY[is] *= -R * R;
      }
      break;
    }
  }

  AMREX_GPU_HOST_DEVICE
  Manifold(const EosParm<Manifold>* eparm)
  {
    eosparm = eparm;

    if (eosparm->manf_data->manmodel == ManifoldModel::TABLE) {
      const TabFuncData* tf_data =
        static_cast<const TabFuncData*>(eosparm->manf_data);
      manfunc = new pele::physics::TabFunc(tf_data);
    } else {
      const NNFuncData* nnf_data =
        static_cast<const NNFuncData*>(eosparm->manf_data);
      manfunc = new pele::physics::NNFunc(nnf_data);
    }
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE Manifold()
  {
    amrex::Error("Manifold EOS must be initialized with table data or neural "
                 "network via EosParm");
  }

  AMREX_GPU_HOST_DEVICE
  ~Manifold() { delete manfunc; }

private:
  pele::physics::ManifoldFunc* manfunc;
  const EosParm<Manifold>* eosparm;
};

} // namespace eos
} // namespace physics
} // namespace pele
#endif
