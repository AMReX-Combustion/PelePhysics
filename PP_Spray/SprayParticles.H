#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include <AMReX_Geometry.H>
#include <AMReX_Particles.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Amr.H>
#include <AMReX_IntVect.H>
#include <AMReX_Gpu.H>
#include "SprayFuelData.H"
#include "EOS.H"

#ifdef USE_SPRAY_SOA
#define NSR_SPR 0
#define NSI_SPR 0
#define NAR_SPR AMREX_SPACEDIM+3+SPRAY_FUEL_NUM
#define NAI_SPR 0
#else
#define NSR_SPR AMREX_SPACEDIM+3+SPRAY_FUEL_NUM
#define NSI_SPR 0
#define NAR_SPR 0
#define NAI_SPR 0
#endif

enum SprayComps {
  heat_tran = 0,
  mass_tran,
  mom_tran,
  pstateVel,
  pstateT,
  pstateRho,
  pstateDia,
  pstateY,
  rhoIndx,
  momIndx,
  engIndx,
  utempIndx,
  specIndx,
  NumIndx
};

class MyParIter
  : public amrex::ParIter<NSR_SPR,NSI_SPR,NAR_SPR,NAI_SPR>
{
public:
  using amrex::ParIter<NSR_SPR,NSI_SPR,NAR_SPR,NAI_SPR>::ParIter;

  std::array<RealVector, NAR_SPR>& GetAttribs() {
    return GetStructOfArrays().GetRealData();
  }

  RealVector& GetAttribs(int comp) {
    return GetStructOfArrays().GetRealData(comp);
  }
};

class MyParConstIter
  : public amrex::ParConstIter<NSR_SPR,NSI_SPR,NAR_SPR,NSI_SPR>
{
public:
  using amrex::ParConstIter<NSR_SPR,NSI_SPR,NAR_SPR,NSI_SPR>::ParConstIter;

  const std::array<RealVector, NAR_SPR>& GetAttribs() const {
    return GetStructOfArrays().GetRealData();
  }

  const RealVector& GetAttribs(int comp) const {
    return GetStructOfArrays().GetRealData(comp);
  }
};

class SprayParticleContainer
  : public amrex::AmrParticleContainer<NSR_SPR,NSI_SPR,NAR_SPR,NAI_SPR>
{
public:

  using PairIndex = std::pair<int, int>;
  using HostVectReal = amrex::Gpu::HostVector<amrex::Real>;
  using HostVectInt = amrex::Gpu::HostVector<int>;
  using GpuArrayInt = amrex::GpuArray<int, SprayComps::NumIndx>;

  SprayParticleContainer(amrex::Amr* amr, amrex::BCRec* _phys_bc, const amrex::Real parcelSize)
    : amrex::AmrParticleContainer<NSR_SPR,NSI_SPR,NAR_SPR,NAI_SPR>(amr),
      m_setFuelData(false),
      m_injectVel(0.),
      m_parcelSize(parcelSize),
      sub_cycle(amr->subCycle()),
      phys_bc(_phys_bc),
      m_sprayRefT(-1.)
  {
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
                     amrex::Real         ref_T,
                     GpuArrayInt&        sprayIndices)
  {
    HostVectReal fuelEnth(NUM_SPECIES);
    EOS::T2Hi(ref_T, fuelEnth.dataPtr());
    m_fuelData.build(critT, boilT, fuelCp, fuelEnth, fuelLatent, fuelIndx);
    m_sprayRefT = ref_T;
    for (int i = 0; i < SprayComps::NumIndx; ++i) {
      AMREX_ASSERT(sprayIndices[i] >= 0);
      m_sprayIndx[i] = sprayIndices[i];
    }
    m_setFuelData = true;
  }

  ~SprayParticleContainer() {}

  ///
  /// Set the value of particle state for all partilces on a level
  ///
  //void SetAll (amrex::Real val, int pstate_idx, int lev);

  ///
  /// Compute a maximum time step based on the particle velocities and a particle CFL number
  ///
  amrex::Real estTimestep(int level, amrex::Real cfl) const;

  ///
  /// Set the number of particles per parcel
  ///
  inline void setParcelSize(const amrex::Real parcelSize)
  {
    m_parcelSize = parcelSize;
  }

  ///
  /// Update the particles but do not move particles
  ///
  void moveKick(amrex::MultiFab& state, amrex::MultiFab& source,
                const int level, const amrex::Real& dt,
                const amrex::Real time, const bool isVirtualPart,
                const bool isGhostPart, const int state_ghosts,
                const int source_ghosts, amrex::MultiFab* u_mac = nullptr);

  ///
  /// Update and move positions of the particles
  ///
  void moveKickDrift(amrex::MultiFab& state, amrex::MultiFab& source,
                     const int level, const amrex::Real& dt,
                     const amrex::Real time, const bool isVirtualPart,
                     const bool isGhostPart, const int state_ghosts,
                     const int source_ghosts, const bool do_move,
                     const int where_width, amrex::MultiFab* u_mac = nullptr);

  ///
  /// Update particles
  ///
  void updateParticles(const int&         level,
                       amrex::MultiFab&   state,
                       amrex::MultiFab&   source,
                       const amrex::Real& flow_dt,
                       const amrex::Real& time,
                       const int          state_ghosts,
                       const int          source_ghosts,
                       const bool         do_move,
                       amrex::MultiFab*   u_mac);

  ///
  /// Sum boundary and copy temporary source term
  ///
  inline void transferSource(const int&       source_ghosts,
                             const int&       level,
                             amrex::MultiFab& tmp_source,
                             amrex::MultiFab& source)
  {
    BL_PROFILE("ParticleContainer::transferSource()");
    const int nghost = amrex::min(source.nGrow(), source_ghosts);
    tmp_source.SumBoundary(Geom(level).periodicity());
    amrex::MultiFab::Add(source, tmp_source, 0, 0, source.nComp(), nghost);
    tmp_source.setVal(0.);
  }

  ///
  /// Insert particles according to some criterion
  ///
  bool insertParticles(amrex::Real time, amrex::Real dt, int nstep, int level, int finest_level);

  bool injectParticles(amrex::Real time, amrex::Real dt, int nstep, int level, int finest_level);

  void InitSprayParticles();

private:
  bool m_setFuelData;
  amrex::Real m_injectVel;
  // The number of spray droplets per computational particle
  amrex::Real m_parcelSize;
  ///
  /// This defines reflect_lo and reflect_hi from phys_bc
  ///
  void init_bcs();

  bool sub_cycle;
  amrex::BCRec* phys_bc;
  bool reflect_lo[AMREX_SPACEDIM];
  bool reflect_hi[AMREX_SPACEDIM];
  SprayDataContainer m_fuelData;
  amrex::Real m_sprayRefT;
  amrex::GpuArray<int, SprayComps::NumIndx> m_sprayIndx;
};

#endif

