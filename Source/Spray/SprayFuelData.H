#ifndef SPRAYFUELDATA_H
#define SPRAYFUELDATA_H

#include "PelePhysics.H"
#include "UnitConversions.H"
#include <AMReX_RealVect.H>

// Spray flags and indices
struct SprayComps
{
  static const int pstateVel = 0; // Particle indices
  static const int pstateT = AMREX_SPACEDIM;
  static const int pstateDia = pstateT + 1;
  static const int pstateY = pstateDia + 1;
  static const int pstateNumDens = pstateY + SPRAY_FUEL_NUM;
  static const int pstateN0 = pstateNumDens + 1;
  static const int pstateBM1 = pstateN0 + 1;  // Breakup model variables
  static const int pstateBM2 = pstateBM1 + 1; // Breakup model variables
  static const int pstateFilmHght = pstateBM2 + 1;
  static const int pstateNum = pstateFilmHght + 1;
  int rhoIndx; // Component indices for conservative variable data structure
  int momIndx;
  int engIndx;
  int utempIndx;
  int specIndx;
  int rhoSrcIndx; // Component indices for spray source data structure
  int momSrcIndx;
  int engSrcIndx;
  int specSrcIndx;
};

enum splash_breakup {
  no_change = 0,
  breakup_KH,
  breakup_RT,
  breakup_TAB,
  splash_dry_splash,
  splash_thermal_breakup,
  splash_wet_splash
};

namespace c2m = pele::physics::utilities::cgs2mks;

// Units for sprays
struct SprayUnits
{
#ifdef PELELM_USE_SPRAY
  // For converting CGS to MKS
  static constexpr amrex::Real ru_conv = c2m::Energy(1.0);
  static constexpr amrex::Real eng_conv = c2m::H(1.0);
  static constexpr amrex::Real rho_conv = c2m::Rho(1.0);
  static constexpr amrex::Real mass_conv = c2m::Mass(1.0);
  static constexpr amrex::Real rhod_conv =
    c2m::Rho(1.0) * c2m::Alpha(1.0); // rho D conversion
  static constexpr amrex::Real mu_conv = c2m::Mu(1.0);
  static constexpr amrex::Real lambda_conv = c2m::Lambda(1.0);
  static constexpr amrex::Real pres_conv = c2m::P(1.0);
  static constexpr amrex::Real len_conv = c2m::Length(1.0);
#else
  static constexpr amrex::Real ru_conv = 1.;
  static constexpr amrex::Real eng_conv = 1.;
  static constexpr amrex::Real rho_conv = 1.;
  static constexpr amrex::Real mass_conv = 1.;
  static constexpr amrex::Real rhod_conv = 1.;
  static constexpr amrex::Real mu_conv = 1.;
  static constexpr amrex::Real lambda_conv = 1.;
  static constexpr amrex::Real pres_conv = 1.;
  static constexpr amrex::Real len_conv = 1.;
#endif
  static constexpr amrex::Real min_mass = 1.E-13 * mass_conv;
};

// Structure that contains the interpolated gas phase state
// and gas phase source values
struct GasPhaseVals
{
  amrex::RealVect vel_fluid;
  amrex::Real T_fluid;
  amrex::Real rho_fluid;
  amrex::Real p_fluid;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_fluid;
  amrex::Real mw_mix;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> mw;
  amrex::RealVect fluid_mom_src;
  amrex::Real fluid_mass_src;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> fluid_Y_dot;
  amrex::Real fluid_eng_src;

  AMREX_GPU_DEVICE AMREX_FORCE_INLINE void reset()
  {
    fluid_mom_src = amrex::RealVect::TheZeroVector();
    vel_fluid = amrex::RealVect::TheZeroVector();
    T_fluid = 0.;
    rho_fluid = 0.;
    fluid_eng_src = 0.;
    fluid_mass_src = 0.;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_fluid[n] = 0.;
    }
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
      fluid_Y_dot[n] = 0.;
    }
  }
};

// Structure containing values for the liquid sprays
struct SprayData
{
  bool mass_trans = true;   // If evaporation is on
  bool mom_trans = true;    // If momentum transfer is on
  bool fixed_parts = false; // If particles are fixed in place
  bool do_splash = false;
  bool liquid_spec_share_gas_dep = false;
  int do_breakup = 0; // 0 - no breakup modeling, 1 - TAB model, 2 - KHRT model
  // Min cell volume fraction to add sources to
  amrex::Real min_eb_vfrac = 0.05;
  amrex::Real ref_T;
  amrex::Real sigma = -1.; // Surface tension
  amrex::Real wall_T = -1.;
  amrex::Real theta_c = -1.; // Contact angle for wall film
  // If particle is updated half dt or whole dt
  amrex::Real dtmod = 0.5;
  amrex::RealVect body_force = amrex::RealVect::TheZeroVector();
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> boilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latent;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> ref_latent;
  // 3 coefficients for Antoine equation and conversion to appropriate units
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> psat_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> rho_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> lambda_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> mu_coef;
  amrex::GpuArray<int, SPRAY_FUEL_NUM> indx = {{-1}};
  amrex::GpuArray<int, SPRAY_FUEL_NUM> dep_indx = {{-1}};

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real rhoL(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = rho_coef[4 * spf];
    amrex::Real b = rho_coef[4 * spf + 1];
    amrex::Real c = rho_coef[4 * spf + 2];
    amrex::Real d = rho_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real lambdaL(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = lambda_coef[4 * spf];
    amrex::Real b = lambda_coef[4 * spf + 1];
    amrex::Real c = lambda_coef[4 * spf + 2];
    amrex::Real d = lambda_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real muL(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = mu_coef[4 * spf];
    amrex::Real b = mu_coef[4 * spf + 1];
    amrex::Real c = mu_coef[4 * spf + 2];
    amrex::Real d = mu_coef[4 * spf + 3];
    return a + ((d / T + c) / T + b) / T;
  }

  // Estimate the boil temperature
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void calcBoilT(const GasPhaseVals& gpv, amrex::Real* cBoilT) const
  {
    SprayUnits SPU;
    amrex::Real RU = pele::physics::Constants::RU * SPU.ru_conv;
    amrex::Real PATM = pele::physics::Constants::PATM * SPU.pres_conv;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fspec = indx[spf];
      const amrex::Real mw_fuel = gpv.mw[fspec];
      // Since we only know the latent heat at the reference temperature,
      // modify Watsons power law to find latent heat at boiling conditions
      amrex::Real Hboil_ref =
        ref_latent[spf] *
        std::pow((critT[spf] - ref_T) / (critT[spf] - boilT[spf]), -0.38);
      // Estimate the boiling temperature at the gas phase pressure using
      // Clasius-Clapeyron relation
      cBoilT[spf] =
        1. / (std::log(PATM / gpv.p_fluid) * RU / (Hboil_ref * mw_fuel) +
              1. / boilT[spf]);
      cBoilT[spf] = amrex::min(critT[spf], cBoilT[spf]);
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real psat(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = psat_coef[4 * spf];
    amrex::Real b = psat_coef[4 * spf + 1];
    amrex::Real c = psat_coef[4 * spf + 2];
    amrex::Real d = psat_coef[4 * spf + 3];
    return d * std::pow(10., a - b / (T + c));
  }
};

#endif
