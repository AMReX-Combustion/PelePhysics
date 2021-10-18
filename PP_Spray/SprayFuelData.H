#ifndef _SPRAYFUELDATA_H_
#define _SPRAYFUELDATA_H_

#include "PelePhysics.H"
#include <AMReX_RealVect.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuMemory.H>

using namespace amrex;

// Spray flags and indices
struct SprayComps
{
  int mass_tran;
  int mom_tran;
  int pstateVel = 0; // Particle indices
  int pstateT = AMREX_SPACEDIM;
  int pstateDia = pstateT + 1;
  int pstateY = pstateDia + 1;
  int pstateVol = pstateVel; // Used when particle is wall film
  int pstateHt = pstateDia;  // Used when particle is wall film
  int rhoIndx; // Component indices for conservative variable data structure
  int momIndx;
  int engIndx;
  int utempIndx;
  int specIndx;
  int rhoSrcIndx; // Component indices for spray source data structure
  int momSrcIndx;
  int engSrcIndx;
  int specSrcIndx;
  int wf_vol = 0; // Wall film indices
  int wf_temp = 1;
  int wf_ht = 2; // Film height
  int wf_Y = 3;
  int wf_num = wf_Y + SPRAY_FUEL_NUM;
};

// Units for sprays
struct SprayUnits
{
#ifdef SPRAY_PELE_LM
  // For converting CGS to MKS
  Real ru_conv = 1.E-7;
  Real eng_conv = 1.E-4;
  Real rho_conv = 1.E3;
  Real mass_conv = 1.E-3;
  Real rhod_conv = 0.1; // rho D conversion
  Real mu_conv = 0.1;
  Real lambda_conv = 1.E-5;
  Real pres_conv = 0.1;
#else
  Real ru_conv = 1.;
  Real eng_conv = 1.;
  Real rho_conv = 1.;
  Real mass_conv = 1.;
  Real rhod_conv = 1.;
  Real mu_conv = 1.;
  Real lambda_conv = 1.;
  Real pres_conv = 1.;
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
  GpuArray<Real, SPRAY_FUEL_NUM> Y_refl;
  RealVect pos_refl;
  RealVect norm;    // Normal to surface
  RealVect tanBeta; // Tangent vector along particle path
  RealVect tanPsi;  // Tangent normal to particle path
  Real omega;       // Used to compute the deviation angle for particle
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
  GpuArray<Real, NUM_SPECIES> Y_fluid;
  GpuArray<Real, NUM_SPECIES> X_fluid;
  Real mw_mix;
  GpuArray<Real, NUM_SPECIES> mw_fluid;
  GpuArray<Real, NUM_SPECIES> invmw;
  RealVect fluid_mom_src;
  Real fluid_mass_src;
  GpuArray<Real, SPRAY_FUEL_NUM> fluid_Y_dot;
  Real fluid_eng_src;

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE GasPhaseVals()
    : vel_fluid(RealVect::TheZeroVector()),
      T_fluid(0.),
      rho_fluid(0.),
      fluid_mom_src(RealVect::TheZeroVector()),
      fluid_mass_src(0.),
      fluid_eng_src(0.)
  {
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_fluid[n] = 0.;
    }
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
      fluid_Y_dot[n] = 0.;
    }
  }

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void define()
  {
    SprayUnits SPU;
    mw_mix = 0.;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      mw_mix += Y_fluid[n] * invmw[n];
    }
    p_fluid =
      rho_fluid * pele::physics::Constants::RU * mw_mix * T_fluid * SPU.ru_conv;
    mw_mix = 1. / mw_mix;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      X_fluid[n] = Y_fluid[n] * invmw[n] * mw_mix;
    }
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
  GpuArray<Real, SPRAY_FUEL_NUM> ref_latent;
  GpuArray<Real, SPRAY_FUEL_NUM> rho;
  GpuArray<Real, SPRAY_FUEL_NUM> mu;
  GpuArray<Real, SPRAY_FUEL_NUM> lambda;
  // 3 coefficients for Antoine equation and conversion to appropriate units
  GpuArray<Real, SPRAY_FUEL_NUM * 4> psat_coef;
  GpuArray<int, SPRAY_FUEL_NUM> indx;

  void build(const SprayData& fdat)
  {
    // Convert input values from CGS to MKS for PeleLM
    num_ppp = fdat.num_ppp;
    ref_T = fdat.ref_T;
    sigma = fdat.sigma;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      critT[spf] = fdat.critT[spf];
      boilT[spf] = fdat.boilT[spf];
      cp[spf] = fdat.cp[spf];
      latent[spf] = fdat.latent[spf];
      ref_latent[spf] = fdat.ref_latent[spf];
      rho[spf] = fdat.rho[spf];
      lambda[spf] = fdat.lambda[spf];
      mu[spf] = fdat.mu[spf];
      for (int cf = 0; cf < 4; ++cf) {
        psat_coef[4 * spf + cf] = fdat.psat_coef[4 * spf + cf];
      }
      indx[spf] = fdat.indx[spf];
    }
  }
};

#endif
