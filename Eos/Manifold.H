#ifndef _MANIFOLD_H_
#define _MANIFOLD_H_

#include <memory>

#ifdef USE_MANIFOLD_TABLE
#include "Table.H"
#endif

#ifdef USE_MANIFOLD_NN
#include "NeuralNet.H"
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
  static void RY2R(const amrex::Real RY[], amrex::Real &R)
  {
    R = RY[NUM_SPECIES-1];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void RY2RRinvY(const amrex::Real RY[], amrex::Real &R, amrex::Real &Rinv, amrex::Real Y[] )
  {
    RY2R(RY,R);
    Rinv = 1.0/R;
    for (int is = 0; is < NUM_SPECIES; is++) {
      Y[is] = RY[is] * Rinv;
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void molecular_weight(amrex::Real* /*mw[]*/)
  {
    amrex::Error("molecular_weight for species does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void inv_molecular_weight(amrex::Real* /*imw[]*/)
  {
    amrex::Error("inv_molecular_weight for species does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(const amrex::Real* /*Y*/, amrex::Real& /*WBAR*/)
  {
    amrex::Error("Y2WBAR is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real* /*Ei[]*/)
  {
    amrex::Error("RTY2Ei does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(const amrex::Real /*T*/, amrex::Real Hi[])
  {
    // TODO: HACK HACK HACK: FIXME (need to remove enthalpy calculations from PeleLM)
    for (int is = 0; is < NUM_SPECIES; is++){
      Hi[is] = 0.0;
    }
    //amrex::Error("T2Hi does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real* /*Hi[]*/)
  {
    amrex::Error("RTY2Hi does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& /*Cv*/)
  {
    amrex::Error("TY2Cv is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cv(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& /*Cv*/)
  {
    amrex::Error("RTY2Cv is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& Cp)
  {
    // TODO: HACK HACK HACK: FIXME (need to remove enthalpy calculations from PeleLM)
    Cp = 100.0;
    //amrex::Error("TY2Cp is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cp(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& /*Cp*/)
  {
    amrex::Error("RTY2Cp is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(const amrex::Real /*R*/, const amrex::Real /*E*/, const amrex::Real* /*Y[]*/, amrex::Real& /*T*/)
  {
    amrex::Error("REY2T is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(const amrex::Real /*H*/, const amrex::Real Y[], amrex::Real& T)
  {
    manfunc->get_value(eosparm->idx_T, Y, T);
    // amrex::Error("HY2T is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(const amrex::Real /*R*/, const amrex::Real /*H*/, const amrex::Real Y[], amrex::Real& T)
  {
    manfunc->get_value(eosparm->idx_T, Y, T);
    // amrex::Error("RHY2T is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real& H)
  {
    H = 0.0;
    // TODO: Eliminate this call entirely?
    // amrex::Error("TY2H is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void X2Y(const amrex::Real X[], amrex::Real Y[])
  {
    // TODO: FIXME HACK HACK HACK
    // just return Y as a hack for PMF, this should really not be needed anywhere
    for (int is = 0; is < NUM_SPECIES; is++){
      Y[is] = X[is];
    }
    //amrex::Error("X2Y does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2X(const amrex::Real* /*Y[]*/, amrex::Real* /*X[]*/)
  {
    amrex::Error("Y2X does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2WDOT(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real Y[], amrex::Real WDOT[])
  {
    // Look up species sources in table if they exist, otherwise they are 0
    for(int is = 0; is < NUM_SPECIES-1; is++)
    {
      if(eosparm->idx_Wdot[is] >= 0) 
      {
	      manfunc->get_value(eosparm->idx_Wdot[is], Y, WDOT[is]);
      }
      else
      {
        manfunc->calculate_Wdot(is, Y, WDOT[is]);
      }
    }

    // Density source is 0
    WDOT[NUM_SPECIES-1] = 0.0;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real* /*JAC[]*/, const int /*HP*/)
  {
    amrex::Error("RTY2JAC is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2C(const amrex::Real /*R*/, const amrex::Real /*T*/, const amrex::Real* /*Y[]*/, amrex::Real* /*JAC[]*/)
  {
    amrex::Error("RTY2C does not have significance for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(const amrex::Real R, const amrex::Real /*T*/, const amrex::Real Y[], amrex::Real& P)
  {
    // Get value of Rho consistent with other state variables
    amrex::Real rho_from_table;
    manfunc->get_value(eosparm->idx_density, Y, rho_from_table);

    // Return pressure is nominal pressure scaled due to density mismatch
    P = eosparm->Pnom_cgs * R/rho_from_table;

    // amrex::Error("RTY2P is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(const amrex::Real /*P*/, const amrex::Real Y[], const amrex::Real /*T*/, amrex::Real& R)
  {
    manfunc->get_value(eosparm->idx_density, Y, R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2dRdY(const amrex::Real Y[], amrex::Real dRdY[])
  {
    manfunc->get_derivs(eosparm->idx_density, Y, dRdY);
    // amrex::Error("Y2dRdY is not yet implemented for Manifold EOS");
  }

  AMREX_GPU_HOST_DEVICE
  Manifold(const EosParm<Manifold>* eparm)
  {
#ifdef USE_MANIFOLD_TABLE
    TabFuncParams::TabFuncData* tf_data =
        static_cast<TabFuncParams::TabFuncData*>(eosparm->manf_data);
    manfunc.reset(new pele::physics::TabFunc(tf_data));
#elif defined(USE_MANIFOLD_NN)
    NNFuncParams::NNFuncData* nnf_data =
        static_cast<NNFuncParams::NNFuncData*>(eosparm->manf_data);
    manfunc.reset(new pele::physics::NNFunc(nnf_data));
#else
    amrex::Error("Invalid option for manifold model.");
#endif
    eosparm = eparm;
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE
  Manifold()
  {
    amrex::Error("Manifold EOS must be initialized with table data or neural network via EosParm");
  }

private:
  std::unique_ptr<pele::physics::ManifoldFunc> manfunc;
  const EosParm<Manifold>* eosparm;
};


}
}
}
#endif
