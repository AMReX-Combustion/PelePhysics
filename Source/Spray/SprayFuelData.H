#ifndef SPRAYFUELDATA_H
#define SPRAYFUELDATA_H

#include "PelePhysics.H"
#include "SprayPropertiesPeleMP.H"
#include "SprayPropertiesGCM.H"
#include <AMReX_RealVect.H>

// Things get tricky for Spray + Manifold because for Manifold EOS, the
// "species" are actually the manifold parameters, but we need the actual
// species for vapor/liquid calculations in the spray module. For manifold, this
// is set up so we only need to lookup gas phase species that correspond to
// lqiuid phase species, rather than all possible gas phase species. We also set
// up macros to define 0 length arrays if variables are unneeded when not using
// Manifold EOS.
#ifdef USE_MANIFOLD_EOS
#define NUM_CHEM_SPECIES SPRAY_FUEL_NUM
#define NUM_CHEM_SPECIES_MANI_ONLY NUM_CHEM_SPECIES
#define NUM_SPECIES_MANI_ONLY NUM_SPECIES
#define SPRAY_FUEL_NUM_MANI_ONLY SPRAY_FUEL_NUM
#else
#define NUM_CHEM_SPECIES NUM_SPECIES
#define NUM_CHEM_SPECIES_MANI_ONLY 0
#define NUM_SPECIES_MANI_ONLY 0
#define SPRAY_FUEL_NUM_MANI_ONLY 0
#endif

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

// Structure that contains the interpolated gas phase state
// and gas phase source values
struct GasPhaseVals
{
  amrex::RealVect vel_fluid;
  amrex::Real T_fluid;
  amrex::Real rho_fluid;
  amrex::Real p_fluid;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_fluid{0.0};
  // for manifold Y_fluid is manifold parameters, Yspec_fluid is chemical
  // species - Yspec_fluid is unused for detailed chem
  amrex::GpuArray<amrex::Real, NUM_CHEM_SPECIES_MANI_ONLY> Yspec_fluid{0.0};
  amrex::Real mw_mix;
  amrex::GpuArray<amrex::Real, NUM_CHEM_SPECIES> mw{0.0};
  amrex::RealVect fluid_mom_src;
  amrex::Real fluid_mass_src;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> fluid_Y_dot{0.0};
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
    for (int n = 0; n < NUM_SPECIES; ++n) {
      fluid_Y_dot[n] = 0.;
    }
  }
};

// Structure containing values for the liquid sprays
struct SprayData
{
  int verbose = 0;          // Same as m_verbose in the SprayParticleContainer
  bool mass_trans = true;   // If evaporation is on
  bool mom_trans = true;    // If momentum transfer is on
  bool fixed_parts = false; // If particles are fixed in place
  bool do_splash = false;
  int do_breakup = 0; // 0 - no breakup modeling, 1 - TAB model, 2 - KHRT model
  // Min cell volume fraction to add sources to
  amrex::Real min_eb_vfrac = 0.05;
  amrex::Real wall_T = -1.;
  amrex::Real theta_c = -1.; // Contact angle for wall film
  // If particle is updated half dt or whole dt
  amrex::Real dtmod = 0.5;
  amrex::RealVect body_force = amrex::RealVect::TheZeroVector();
  amrex::GpuArray<int, SPRAY_FUEL_NUM> dep_indx{};
  amrex::GpuArray<int, SPRAY_FUEL_NUM> L_col{};
  amrex::GpuArray<int, NUM_CHEM_SPECIES + 1> L_row{};
  int N_pc = 0;
  // We will only use the first N_pc entries, but don't know that number a
  // priori so allocate NUM_SPECIES which is always >= N_pc
  amrex::GpuArray<int, NUM_CHEM_SPECIES> pc_indx;

  // Arrays only used for Manifold EOS
  amrex::GpuArray<int, SPRAY_FUEL_NUM_MANI_ONLY> dep_manifold_indx{};
  amrex::GpuArray<int, NUM_CHEM_SPECIES_MANI_ONLY> pc_manifold_indx{};

  // Instance of the liquid properties for the specified model
  pele::physics::SprayProps::LiqPropType liqprops;
  // adding a eosparm here to be used with Spray and manifold
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparm;

  // Estimate the boil temperature
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void calcBoilT(const GasPhaseVals& gpv, amrex::Real* cBoilT) const
  {
    const amrex::Real RU = pele::physics::Constants::RU * SprayUnits::ru_conv;
    const amrex::Real PATM =
      pele::physics::Constants::PATM * SprayUnits::pres_conv;

    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const amrex::Real mw_li = liqprops.mw_i(spf);

      // Use correlation to avoid calling eos.T2Hi for each boilT_i(spf)
      const amrex::Real Hboil_ref =
        liqprops.latent_correlation_i(liqprops.boilT_i(spf), spf);

      // Estimate the boiling temperature at the gas phase pressure using
      // Clasius-Clapeyron relation
      cBoilT[spf] =
        1. / (std::log(PATM / gpv.p_fluid) * RU / (Hboil_ref * mw_li) +
              1. / liqprops.boilT_i(spf));
      cBoilT[spf] = amrex::min(liqprops.critT_i(spf), cBoilT[spf]);
    }
  }

  // Compute the dot product of the i-th row of a binary CSR matrix L with
  // vector x.
  template <
    typename IntArrayLike,
    typename IntArrayLike2,
    typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real binary_csr_rowdotprod(
    int i,
    const IntArrayLike& row_ptr,
    const IntArrayLike2& col_idx,
    const RealArrayLike& x) const
  {
    amrex::Real sum = 0.;
    for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
      sum += x[col_idx[j]];
    }
    return sum;
  }

  // Return an ordered vector of unique nonzero row indices in a CSR matrix.
  template <typename IntArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
  binary_csr_nonzerorows(const IntArrayLike& row_ptr, int* nonzero_rows)
  {
    int num_rows = row_ptr.size() - 1;
    int nnzrows = 0;
    for (int i = 0; i < num_rows; ++i) {
      if (row_ptr[i + 1] > row_ptr[i]) {
        nonzero_rows[nnzrows++] = i;
      }
    }
    return nnzrows;
  }
};

#endif
