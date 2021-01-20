#ifndef _SPRAYFUELDATA_H_
#define _SPRAYFUELDATA_H_

#include <AMReX_Gpu.H>
#include "EOS.H"

using namespace amrex;

// Spray flags and indices
struct SprayComps {
  int heat_tran;
  int mass_tran;
  int mom_tran;
  int pstateVel; // Particle indices
  int pstateT;
  int pstateDia;
  int pstateY;
  int rhoIndx; // Gas phase indices
  int momIndx;
  int engIndx;
  int utempIndx;
  int specIndx;
};

// Units for sprays
struct SprayUnits {
#ifdef SPRAY_PELE_LM
  Real vel_conv = 100.;  // Turn m/s to cm/s
  Real rho_conv = 0.001; // Turn kg/m^3 to g/cm^3
  Real pos_conv = 0.01;  // Turn cm to m for updating position
  Real eng_conv = 1.E4; // For converting enthalpy to CGS
  // This makes no sense, conversions should be independent
  // of dimensions but numerical tests show this isn't the case
#if AMREX_SPACEDIM == 2
  Real mom_src_conv = 1.E-3;
  Real mass_src_conv = 1.E-1;
  Real eng_src_conv = 1.E-5;
#elif AMREX_SPACEDIM == 3
  Real mom_src_conv = 1.E-5;
  Real mass_src_conv = 1.E-3;
  Real eng_src_conv = 1.E-7;
#endif
#else
  Real vel_conv = 1.;
  Real pos_conv = 1.;
  Real rho_conv = 1.;
  Real eng_conv = 1.;
  Real mom_src_conv = 1.;
  Real mass_src_conv = 1.;
  Real eng_src_conv = 1.;
#endif
};

// Structure that contains the values for secondary
// droplets that form during impingement
struct SprayRefl
{
  Real Unorm;
  Real dt_pp; // Time remaining to convect reflected drops
  Real dia_refl;
  int Ns_refl = 0;
  GpuArray<Real,SPRAY_FUEL_NUM> Y_refl;
  RealVect norm; // Normal to surface
  RealVect tanBeta; // Tangent vector along particle path
  RealVect tanPsi; // Tangent normal to particle path
  Real omega; // Used to compute the deviation angle for particle
  Real expomega;
  Real beta_mean; // Mean ejection angle in degrees
  Real beta_stdv; // Standard dev ejection angle in degrees
};

// Structure that contains the interpolated gas phase state
// and gas phase source values
struct GasPhaseVals
{
  RealVect vel_fluid;
  Real T_fluid;
  Real rho_fluid;
  Real p_fluid;
  GpuArray<Real,NUM_SPECIES> Y_fluid;
  Real mw_mix;
  GpuArray<Real,NUM_SPECIES> mw_fluid;
  GpuArray<Real,NUM_SPECIES> invmw;
  RealVect fluid_mom_src;
  Real fluid_mass_src;
  GpuArray<Real,SPRAY_FUEL_NUM> fluid_Y_dot;
  Real fluid_eng_src;

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  GasPhaseVals (const RealVect in_vel,
                const Real     in_T,
                const Real     in_rho,
                const Real*    in_Y,
                const Real*    in_mw,
                const Real*    in_invmw)
    : vel_fluid(in_vel),
      T_fluid(in_T),
      rho_fluid(in_rho),
      fluid_mom_src(RealVect::TheZeroVector()),
      fluid_eng_src(0.),
      fluid_mass_src(0.)
  {
    mw_mix = 0.;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_fluid[n] = in_Y[n];
      mw_fluid[n] = in_mw[n];
      invmw[n] = in_invmw[n];
      mw_mix += Y_fluid[n]*invmw[n];
    }
    p_fluid = rho_fluid*EOS::RU*mw_mix*T_fluid;
    mw_mix = 1./mw_mix;
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n)
      fluid_Y_dot[n] = 0.;
  }
};

// Structure containing values for the liquid sprays
struct SprayData
{
  Real num_ppp; // Parcel size
  Real ref_T;
  Real sigma; // Surface tension

  SprayData (const Gpu::DeviceVector<Real>& critT,
             const Gpu::DeviceVector<Real>& boilT,
             const Gpu::DeviceVector<Real>& fuelCp,
             const Gpu::DeviceVector<Real>& fuelLatent,
             const Gpu::DeviceVector<Real>& fuelRho,
             const Gpu::DeviceVector<Real>& fuelMu,
             const Gpu::DeviceVector<int>&  fuelIndx,
             const Real in_num_ppp,
             const Real in_ref_T,
             const Real in_sigma)
    : m_critT_ptr(critT.dataPtr()),
      m_boilT_ptr(boilT.dataPtr()),
      m_fuelCp_ptr(fuelCp.dataPtr()),
      m_fuelLatent_ptr(fuelLatent.dataPtr()),
      m_fuelRho_ptr(fuelRho.dataPtr()),
      m_fuelMu_ptr(fuelMu.dataPtr()),
      m_fuelIndx_ptr(fuelIndx.dataPtr()),
      num_ppp(in_num_ppp),
      ref_T(in_ref_T),
      sigma(in_sigma)
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real critT (const int fuelSpecIndx) const
  {
    return m_critT_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real boilT (const int fuelSpecIndx) const
  {
    return m_boilT_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real cp (const int fuelSpecIndx) const
  {
    return m_fuelCp_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real latent (const int fuelSpecIndx) const
  {
    return m_fuelLatent_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real rho (const int fuelSpecIndx) const
  {
    return m_fuelRho_ptr[fuelSpecIndx];
  }

  // Surface tension of liquid fuel for splash model
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real mu (const int fuelSpecIndx) const
  {
    return m_fuelMu_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int indx (const int fuelSpecIndx) const
  {
    return m_fuelIndx_ptr[fuelSpecIndx];
  }

private:

  const Real * m_critT_ptr;
  const Real * m_boilT_ptr;
  const Real * m_fuelCp_ptr;
  const Real * m_fuelLatent_ptr;
  const Real * m_fuelRho_ptr;
  const Real * m_fuelMu_ptr;
  const int * m_fuelIndx_ptr;
};

class SprayDataContainer
{
public:

  void build (const Gpu::HostVector<Real>& critT,
              const Gpu::HostVector<Real>& boilT,
              const Gpu::HostVector<Real>& fuelCp,
              const Gpu::HostVector<Real>& refFuelH,
              const Gpu::HostVector<Real>& fuelLatent,
              const Gpu::HostVector<Real>& fuelRho,
              const Gpu::HostVector<Real>& fuelMu,
              const Gpu::HostVector<int>&  fuelIndx,
              Real in_num_ppp,
              Real in_ref_T,
              Real in_sigma)
  {
    num_ppp = in_num_ppp;
    ref_T = in_ref_T;
    fuelSigma = in_sigma;
    const int nfspec = SPRAY_FUEL_NUM;
    m_critT.resize(nfspec);
    m_boilT.resize(nfspec);
    m_fuelCp.resize(nfspec);
    m_fuelLatent.resize(nfspec);
    m_fuelRho.resize(nfspec);
    m_fuelMu.resize(nfspec);
    m_fuelIndx.resize(nfspec);
    Gpu::HostVector<Real> tmpFuelLatent(nfspec);
    for (int i = 0; i != nfspec; ++i) {
      const int specIndx = fuelIndx[i];
      tmpFuelLatent[i] = fuelLatent[i] - refFuelH[specIndx];
    }
    Gpu::copy(Gpu::hostToDevice, critT.begin(), critT.end(),
              m_critT.begin());
    Gpu::copy(Gpu::hostToDevice, boilT.begin(), boilT.end(),
              m_boilT.begin());
    Gpu::copy(Gpu::hostToDevice, fuelCp.begin(), fuelCp.end(),
              m_fuelCp.begin());
    Gpu::copy(Gpu::hostToDevice, tmpFuelLatent.begin(), tmpFuelLatent.end(),
              m_fuelLatent.begin());
    Gpu::copy(Gpu::hostToDevice, fuelRho.begin(), fuelRho.end(),
              m_fuelRho.begin());
    Gpu::copy(Gpu::hostToDevice, fuelMu.begin(), fuelMu.end(),
              m_fuelMu.begin());
    Gpu::copy(Gpu::hostToDevice, fuelIndx.begin(), fuelIndx.end(),
              m_fuelIndx.begin());
  }

  SprayData getSprayData()
  {
    return SprayData(m_critT, m_boilT, m_fuelCp, m_fuelLatent, m_fuelRho,
                     m_fuelMu, m_fuelIndx, num_ppp, ref_T, fuelSigma);
  }

private:

  // Vectors with constant reference values
  Gpu::DeviceVector<Real> m_critT;
  Gpu::DeviceVector<Real> m_boilT;
  Gpu::DeviceVector<Real> m_fuelCp;
  Gpu::DeviceVector<Real> m_fuelLatent;
  Gpu::DeviceVector<Real> m_fuelRho;
  Gpu::DeviceVector<Real> m_fuelMu;
  Gpu::DeviceVector<int>  m_fuelIndx;
  Real num_ppp;
  Real ref_T;
  Real fuelSigma; // Surface tension
};

#endif
