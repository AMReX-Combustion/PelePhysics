#ifndef SPRAYFUELDATA_H
#define SPRAYFUELDATA_H

#include "PelePhysics.H"
#include <AMReX_RealVect.H>

// Spray flags and indices
struct SprayComps
{
  int pstateVel = 0; // Particle indices
  int pstateT = AMREX_SPACEDIM;
  int pstateDia = pstateT + 1;
  int pstateY = pstateDia + 1;
  int pstatePb = pstateY + SPRAY_FUEL_NUM;
  int psttePbdot = pstatePb + 1;
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
  breakup,
  splash_splash,
  splash_thermal_breakup
};

// Units for sprays
struct SprayUnits
{
#ifdef PELELM_USE_SPRAY
  // For converting CGS to MKS
  amrex::Real ru_conv = 1.E-7;
  amrex::Real eng_conv = 1.E-4;
  amrex::Real rho_conv = 1.E3;
  amrex::Real mass_conv = 1.E-3;
  amrex::Real rhod_conv = 0.1; // rho D conversion
  amrex::Real mu_conv = 0.1;
  amrex::Real lambda_conv = 1.E-5;
  amrex::Real pres_conv = 0.1;
  amrex::Real len_conv = 0.01;
#else
  amrex::Real ru_conv = 1.;
  amrex::Real eng_conv = 1.;
  amrex::Real rho_conv = 1.;
  amrex::Real mass_conv = 1.;
  amrex::Real rhod_conv = 1.;
  amrex::Real mu_conv = 1.;
  amrex::Real lambda_conv = 1.;
  amrex::Real pres_conv = 1.;
  amrex::Real len_conv = 1.;
#endif
  amrex::Real min_mass = 1.E-13 * mass_conv;
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
  amrex::GpuArray<amrex::Real, NUM_SPECIES> X_fluid;
  amrex::Real mw_mix;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> mw_fluid;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> invmw;
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
  bool mass_trans = true;   // If evaporation is on
  bool mom_trans = true;    // If momentum transfer is on
  bool fixed_parts = false; // If particles are fixed in place
  // Min cell volume fraction to add sources to
  amrex::Real min_eb_vfrac = 0.05;
  amrex::Real num_ppp; // Parcel size
  amrex::Real ref_T;
  amrex::Real sigma = -1.; // Surface tension
  amrex::Real wall_T = -1.;
  // If particle is updated half dt or whole dt
  amrex::Real dtmod = 0.5;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> boilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latent;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> ref_latent;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> rho;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mu;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> lambda;
  // 3 coefficients for Antoine equation and conversion to appropriate units
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> psat_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> rho_coef;
  amrex::GpuArray<int, SPRAY_FUEL_NUM> indx = {{-1}};
  amrex::GpuArray<int, SPRAY_FUEL_NUM> dep_indx = {{-1}};

  void build(const SprayData& fdat)
  {
    mass_trans = fdat.mass_trans;
    mom_trans = fdat.mom_trans;
    fixed_parts = fdat.fixed_parts;
    // Convert input values from CGS to MKS for PeleLM
    num_ppp = fdat.num_ppp;
    ref_T = fdat.ref_T;
    sigma = fdat.sigma;
    wall_T = fdat.wall_T;
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
        rho_coef[4 * spf + cf] = fdat.rho_coef[4 * spf + cf];
      }
      indx[spf] = fdat.indx[spf];
      dep_indx[spf] = fdat.dep_indx[spf];
    }
    dtmod = fdat.dtmod;
    min_eb_vfrac = fdat.min_eb_vfrac;
    T_wall = fdat.T_wall;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real rhoL(const amrex::Real& T, const int spf) const
  {
    amrex::Real rhoL = rho[spf];
    if (rho_coef[4 * spf] != 0.) {
      amrex::Real a = rho_coef[4 * spf];
      amrex::Real b = rho_coef[4 * spf + 1];
      amrex::Real c = rho_coef[4 * spf + 2];
      amrex::Real d = rho_coef[4 * spf + 3];
      rhoL = a + T * (b + T * (c + T * d));
    }
    return rhoL;
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
      const amrex::Real mw_fuel = gpv.mw_fluid[fspec];
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

// Here are different typical jet droplet diameter distributions

// Normal or Gaussian distribution
class NormDist
{
  amrex::Real mean_dia;
  amrex::Real std_dia;

public:
  explicit NormDist(amrex::Real a_mean, amrex::Real a_std)
    : mean_dia(a_mean), std_dia(a_std)
  {
  }
  amrex::Real get_dia() const { return amrex::RandomNormal(mean_dia, std_dia); }
  amrex::Real get_avg_dia() const { return mean_dia; }
};

// Log-normal distribution
class LogNormDist
{
  amrex::Real log_mean;
  amrex::Real log_std;

public:
  ~LogNormDist() = default;
  explicit LogNormDist(amrex::Real a_mean, amrex::Real a_std)
  {
    amrex::Real stdsq = a_std * a_std;
    amrex::Real meansq = a_mean * a_mean;
    log_mean = 2. * std::log(a_mean) - 0.5 * std::log(stdsq + meansq);
    log_std = std::sqrt(
      amrex::max(-2. * std::log(a_mean) + std::log(stdsq + meansq), 0.));
  }
  amrex::Real get_dia() const
  {
    return std::exp(amrex::RandomNormal(log_mean, log_std));
  }
  amrex::Real get_avg_dia() const { return std::exp(log_mean); }
};

// Weibull or Rosin-Rammler distribution
class WeibullDist
{
  amrex::Real mean_dia;
  amrex::Real k;

public:
  ~WeibullDist() = default;
  explicit WeibullDist(amrex::Real a_mean, amrex::Real a_k)
    : mean_dia(a_mean), k(a_k)
  {
  }
  amrex::Real get_dia() const
  {
    amrex::Real fact = -std::log(1. - amrex::Random());
    return mean_dia * std::pow(fact, 1. / k);
  }
  amrex::Real get_avg_dia() const { return mean_dia; }
};

#endif
