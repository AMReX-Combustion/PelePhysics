#ifndef _SPRAYFUELDATA_H_
#define _SPRAYFUELDATA_H_

#include <AMReX_Gpu.H>
#include "EOS.H"

using namespace amrex;

// Spray flags and indices
struct SprayComps
{
  int heat_tran;
  int mass_tran;
  int mom_tran;
  int pstateVel = 0; // Particle indices
  int pstateT = AMREX_SPACEDIM;
  int pstateDia = pstateT + 1;
  int pstateY = pstateDia + 1;
  int rhoIndx; // Gas phase indices
  int momIndx;
  int engIndx;
  int utempIndx;
  int specIndx;
  int wf_vol = 0; // Wall film indices
  int wf_temp = 1;
  int wf_Y = 2;
  int wf_num = wf_Y + SPRAY_FUEL_NUM;
};

// Units for sprays
struct SprayUnits
{
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
  RealVect pos_refl;
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
  GpuArray<Real, SPRAY_FUEL_NUM> critT;
  GpuArray<Real, SPRAY_FUEL_NUM> boilT;
  GpuArray<Real, SPRAY_FUEL_NUM> cp;
  GpuArray<Real, SPRAY_FUEL_NUM> latent;
  GpuArray<Real, SPRAY_FUEL_NUM> rho;
  GpuArray<Real, SPRAY_FUEL_NUM> mu;
  GpuArray<int, SPRAY_FUEL_NUM> indx;

  SprayData (const Real* in_critT,
             const Real* in_boilT,
             const Real* in_fuelCp,
             const Real* in_fuelLatent,
             const Real* in_fuelRho,
             const Real* in_fuelMu,
             const int* in_fuelIndx,
             const Real in_num_ppp,
             const Real in_ref_T,
             const Real in_sigma)
    : num_ppp(in_num_ppp),
      ref_T(in_ref_T),
      sigma(in_sigma)
  {
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      critT[spf] = in_critT[spf];
      boilT[spf] = in_boilT[spf];
      cp[spf] = in_fuelCp[spf];
      latent[spf] = in_fuelLatent[spf];
      rho[spf] = in_fuelRho[spf];
      mu[spf] = in_fuelMu[spf];
      indx[spf] = in_fuelIndx[spf];
    }
  }
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
    for (int i = 0; i != nfspec; ++i) {
      const int specIndx = fuelIndx[i];
      m_fuelLatent[i] = fuelLatent[i] - refFuelH[specIndx];
      m_critT[i] = critT[i];
      m_boilT[i] = boilT[i];
      m_fuelCp[i] = fuelCp[i];
      m_fuelRho[i] = fuelRho[i];
      m_fuelMu[i] = fuelMu[i];
      m_fuelIndx[i] = fuelIndx[i];
    }
  }

  SprayData getSprayData()
  {
    return SprayData(m_critT.data(), m_boilT.data(), m_fuelCp.data(),
                     m_fuelLatent.data(), m_fuelRho.data(), m_fuelMu.data(),
                     m_fuelIndx.data(), num_ppp, ref_T, fuelSigma);
  }

private:

  // Vectors with constant reference values
  Vector<Real> m_critT;
  Vector<Real> m_boilT;
  Vector<Real> m_fuelCp;
  Vector<Real> m_fuelLatent;
  Vector<Real> m_fuelRho;
  Vector<Real> m_fuelMu;
  Vector<int>  m_fuelIndx;
  Real num_ppp;
  Real ref_T;
  Real fuelSigma; // Surface tension
};

#endif
