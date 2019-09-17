#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include <AMReX_Geometry.H>
#include "AMReX_Particles.H"
#include "AMReX_AmrParticles.H"
#include "AMReX_Amr.H"

/*
  Setting the number of components in the particle structs, we must be sure that 
  this is consistent with the definition of the particle type in f90.  We make space
  for position (D), velocity (D), temp (1), diameter (1) and density (1) [here, D is 
  the space dimension].  The library already accounts for the position, so we 
  declear D+3.

  From particle_mod.f90:

  type, bind(C)  :: particle_t
     real(c_real)    :: pos(amrex_spacedim)  !< Position
     real(c_real)    :: vel(amrex_spacedim)  !< Particle velocity
     real(c_real)    :: temp                 !< Particle temperature
     real(c_real)    :: diam                 !< Particle diameter
     real(c_real)    :: density              !< Particle density
     integer(c_int)  :: id
     integer(c_int)  :: cpu
  end type particle_t

 */

#define NSR_SPR AMREX_SPACEDIM+3
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

   SprayParticleContainer(amrex::Amr* amr, amrex::BCRec* _phys_bc)
        : amrex::AmrParticleContainer<NSR_SPR,0>(amr),
        sub_cycle(amr->subCycle()),
        phys_bc(_phys_bc)
	{
            nstate = NSR_SPR;
            init_bcs();
        }

    ~SprayParticleContainer () {}

    void AssignDensity (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, int lev_min = 0, int ncomp = 1,
                        int finest_level = -1)
    {
        //HK: Call to SetAllowParticlesNearBoundary is no longer needed, in fact method no longer exists
        //amrex::AmrParticleContainer<NSR_SPR,0>::SetAllowParticlesNearBoundary (true);
        amrex::AmrParticleContainer<NSR_SPR,0>::AssignDensity(0, mf, lev_min, ncomp, finest_level);
    }

    ///
    /// Set the value of particle state for all partilces on a level
    ///
    void SetAll (amrex::Real val, int pstate_idx, int lev);

    ///
    /// Compute a maximum time step based on the particle velocities and a particle CFL number
    ///
    amrex::Real estTimestep (int level, amrex::Real cfl) const;

    ///
    /// Update the velocity and positions of the particles
    ///
    void moveKickDrift (amrex::MultiFab& state, amrex::MultiFab& source,
                        int level, amrex::Real dt, int tmp_src_width, int where_width);

    ///
    /// Update only the velocity of the particles
    ///
    void moveKick      (amrex::MultiFab& state, amrex::MultiFab& source,
                        int level, amrex::Real dt, int tmp_src_width);

    ///
    /// Insert particles according to some criterion
    ///
    void insertParticles (amrex::Real time, int nstep, int lev);

    void injectParticles (amrex::Real time, int nstep, int lev);

    void InitParticlesUniform(const int& lev, const int& num_ppc);

private:
    int nstate;
    ///
    /// This defines reflect_lo and reflect_hi from phys_bc
    ///
    void init_bcs();

    bool sub_cycle;
    amrex::BCRec* phys_bc;
    int reflect_lo[AMREX_SPACEDIM];
    int reflect_hi[AMREX_SPACEDIM];
};

#endif

