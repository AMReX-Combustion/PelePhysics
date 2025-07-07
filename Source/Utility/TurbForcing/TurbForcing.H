#ifndef TURBFORCING_H
#define TURBFORCING_H

#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <MersenneTwister.H>

namespace pele::physics::turbforcing {

struct TurbForcingParm
{
  int m_verbose = 0;

  // fine tune force scale
  amrex::Real m_force_scale_fudge = 1.0;
  // time offset for forcing function
  amrex::Real m_time_offset = 0.0;
  // target urms for the forcing
  amrex::Real m_urms = 1.0;
  // allow peridodic reproduction in z
  // can be used as a flag i.e. =1 for a factor 2
  // or as the factor itself
  int m_hack_lz = 0;

  // Tuned by Andrew Aspden, change at own risk!
  int m_ff_factor = 4;
  int m_nmodes = 4;
  amrex::Real m_forcing_epsilon = 0.1; // suppress symmetry breaking modes
  int m_spectrum_type = 2;
  int m_moderate_zero_modes = 1;
  int m_mode_start = 0;
  static constexpr int m_array_size = 33;
  static constexpr int m_num_fdarray = 17;

  // Andrew Aspden's reference values for forcing scales
  static constexpr amrex::Real m_flr = 2.64e-2;
  static constexpr amrex::Real m_fvr = 3.27;
  static constexpr amrex::Real m_fts_min = 1.21e-3;
  static constexpr amrex::Real m_fts_max = 2.42e-3;

  amrex::Real m_fsr;
  // computed values for forcing
  amrex::Real m_force_scale;
  amrex::Real m_forcing_time_scale_min;
  amrex::Real m_forcing_time_scale_max;
  amrex::Real m_Lmin;
  amrex::Real m_kappaMax;
  int m_nxmodes;
  int m_nymodes;
  int m_nzmodes;
  amrex::Real m_freqMin;
  amrex::Real m_freqMax;
  amrex::Real m_freqDiff;

  // forcedata will contain num_fdarray arrays of size
  // (0,0,0)(array_size-1,array_size-1,array_size-1)
  amrex::Real* m_forcedata;
};

struct TurbForcing
{
public:
  TurbForcing() = default;

  ~TurbForcing() = default;

  // initialises the parameters
  void init(amrex::GeometryData const& geomdata);

  // builds the forcing term
  // note: a_rho_incompressible is unused if a_incompressible
  // == false.
  void addTurbVelForces(
    amrex::GeometryData const& geomdata,
    const amrex::Box& bx,
    const amrex::Real& a_time,
    amrex::Array4<amrex::Real> const& force,
    amrex::Array4<const amrex::Real> const& rho,
    const int a_incompressible,
    const amrex::Real a_rho_incompressible = -1.0);

private:
  TurbForcingParm m_tfp;
  bool m_turbforcing_initialized = false;
};

} // namespace pele::physics::turbforcing

#endif
