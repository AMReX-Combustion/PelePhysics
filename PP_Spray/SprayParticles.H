#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include <AMReX_Geometry.H>
#include "AMReX_Particles.H"
#include "AMReX_AmrParticles.H"
#include "AMReX_Amr.H"
#include <AMReX_IntVect.H>
#include <AMReX_Gpu.H>
#include "SprayFuelData.H"
#include "EOS.H"

// Define the number of fuel species used in the droplets, if not defined in GNUmakefile
#ifndef SPRAY_FUEL_NUM
#define SPRAY_FUEL_NUM 1
#endif

#define NSR_SPR AMREX_SPACEDIM+3+SPRAY_FUEL_NUM
#define NAR_SPR 0
#define NSI_SPR 0
#define NAI_SPR 0

class SprayParticleContainer
  : public amrex::AmrParticleContainer<NSR_SPR,NSI_SPR,NAR_SPR,NAI_SPR>
{
public:

  using MyParIter = amrex::ParIter<NSR_SPR,0>;
  using MyParConstIter = amrex::ParConstIter<NSR_SPR,0>;
  using PairIndex = std::pair<int, int>;
  using HostVectReal = amrex::Gpu::HostVector<amrex::Real>;
  using HostVectInt = amrex::Gpu::HostVector<int>;

  SprayParticleContainer(amrex::Amr* amr, amrex::BCRec* _phys_bc)
    : amrex::AmrParticleContainer<NSR_SPR,0>(amr),
    sub_cycle(amr->subCycle()),
    phys_bc(_phys_bc)
  {
    nstate = NSR_SPR;
    init_bcs();
  }

  ///
  /// Set fuel data
  ///
  void buildFuelData(const HostVectReal& critT,
		     const HostVectReal& boilT,
		     const HostVectReal& fuelCp,
		     const HostVectReal& fuelLatent,
		     const HostVectInt&  fuelIndx,
		     amrex::Real         ref_T)
  {
    HostVectReal fuelEnth(NUM_SPECIES);
#ifdef PELEC_EOS_FUEGO
    EOS::get_hi(nullptr, ref_T, fuelEnth.dataPtr());
#else
    amrex::Real mass[NUM_SPECIES] = {28.97};
    EOS::get_hi(mass, ref_T, fuelEnth.dataPtr());
#endif
    m_fuelData.build(critT, boilT, fuelCp, fuelEnth, fuelLatent, fuelIndx);
  }

  ~SprayParticleContainer () {}

  ///
  /// Set the value of particle state for all partilces on a level
  ///
  //void SetAll (amrex::Real val, int pstate_idx, int lev);

  ///
  /// Compute a maximum time step based on the particle velocities and a particle CFL number
  ///
  amrex::Real estTimestep (int level, amrex::Real cfl) const;

  ///
  /// Update the particles but do not move particles
  ///
  void moveKick (amrex::MultiFab& state, amrex::MultiFab& source,
		 const int level, const amrex::Real& dt,
		 const amrex::Real time, const int tmp_src_width);

  ///
  /// Update and move positions of the particles
  ///
  void moveKickDrift (amrex::MultiFab& state, amrex::MultiFab& source,
		      const int level, const amrex::Real& dt,
		      const amrex::Real time, const int tmp_src_width,
		      const bool do_move, const int where_width);

  ///
  /// Update particles
  ///
  void updateParticles (const int&         lev,
			amrex::MultiFab&   state,
			amrex::MultiFab&   source,
			const amrex::Real& flow_dt,
			const amrex::Real& time,
			const bool         do_move);

  ///
  /// Insert particles according to some criterion
  ///
  bool insertParticles (amrex::Real time, int nstep, int lev);

  bool injectParticles (amrex::Real time, int nstep, int lev);

  void InitParticlesUniform (amrex::AmrLevel* pelec, const int& lev, const int& num_ppc);

private:
  int nstate;
  ///
  /// This defines reflect_lo and reflect_hi from phys_bc
  ///
  void init_bcs();

  bool sub_cycle;
  amrex::BCRec* phys_bc;
  amrex::IntVect reflect_lo;
  amrex::IntVect reflect_hi;
  SprayFuelData m_fuelData;
};

#endif

