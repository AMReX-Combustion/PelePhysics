#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include "EOS.H"
#include "SprayFuelData.H"
#include <AMReX_Amr.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_Particles.H>
#include <memory>

#ifdef SPRAY_PELE_LM
#include "pelelm_prob.H"
#include "pelelm_prob_parm.H"
#else
#include "prob.H"
#include "prob_parm.H"
#endif

#ifdef USE_SPRAY_SOA
#define NSR_SPR 0
#define NSI_SPR 0
#define NAR_SPR AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM
#define NAI_SPR 0
#else
#define NSR_SPR AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM
#define NSI_SPR 0
#define NAR_SPR 0
#define NAI_SPR 0
#endif

class MyParIter : public amrex::ParIter<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>
{
public:
  using amrex::ParIter<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>::ParIter;
#ifdef USE_SPRAY_SOA
  std::array<RealVector, NAR_SPR>& GetAttribs()
  {
    return GetStructOfArrays().GetRealData();
  }

  RealVector& GetAttribs(int comp)
  {
    return GetStructOfArrays().GetRealData(comp);
  }
#endif
};

class MyParConstIter
  : public amrex::ParConstIter<NSR_SPR, NSI_SPR, NAR_SPR, NSI_SPR>
{
public:
  using amrex::ParConstIter<NSR_SPR, NSI_SPR, NAR_SPR, NSI_SPR>::ParConstIter;
#ifdef USE_SPRAY_SOA
  const std::array<RealVector, NAR_SPR>& GetAttribs() const
  {
    return GetStructOfArrays().GetRealData();
  }

  const RealVector& GetAttribs(int comp) const
  {
    return GetStructOfArrays().GetRealData(comp);
  }
#endif
};

class SprayParticleContainer
  : public amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>
{
public:
  using PairIndex = std::pair<int, int>;
  using HostVectReal = amrex::Gpu::HostVector<amrex::Real>;
  using HostVectInt = amrex::Gpu::HostVector<int>;

  SprayParticleContainer(
    amrex::Amr* amr,
    amrex::BCRec* _phys_bc,
    SprayData fdat,
    SprayComps SPI,
    const amrex::Real parcelSize,
    const amrex::Real wallT)
    : amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>(amr),
      m_injectVel(0.),
      m_parcelSize(parcelSize),
      m_wallT(wallT),
      phys_bc(_phys_bc),
      m_sprayData(nullptr),
      d_sprayData(nullptr),
      m_sprayIndx(SPI)
  {
    m_sprayData = new SprayData{};
    d_sprayData =
      static_cast<SprayData*>(amrex::The_Arena()->alloc(sizeof(SprayData)));
    m_sprayData->build(fdat);
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, m_sprayData, m_sprayData + 1, d_sprayData);
    init_bcs();
  }

  ~SprayParticleContainer()
  {
    delete m_sprayData;
    amrex::The_Arena()->free(d_sprayData);
  }

  void SprayParticleIO(
    const int level,
    const bool is_checkpoint,
    const int write_ascii,
    const std::string& dir,
    const amrex::Vector<std::string>& sprayFuelNames)
  {
    const int pstateVel = m_sprayIndx.pstateVel;
    const int pstateT = m_sprayIndx.pstateT;
    const int pstateDia = m_sprayIndx.pstateDia;
    const int pstateY = m_sprayIndx.pstateY;
    amrex::Vector<std::string> real_comp_names(pstateY + SPRAY_FUEL_NUM);
    AMREX_D_TERM(real_comp_names[pstateVel] = "xvel";
                 , real_comp_names[pstateVel + 1] = "yvel";
                 , real_comp_names[pstateVel + 2] = "zvel";);
    real_comp_names[pstateT] = "temperature";
    real_comp_names[pstateDia] = "diam";
    for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp) {
      real_comp_names[pstateY + sp] = "spray_mf_" + sprayFuelNames[sp];
    }
    amrex::Vector<std::string> int_comp_names;
    Checkpoint(
      dir, "particles", is_checkpoint, real_comp_names, int_comp_names);
    // Here we write ascii information every time we write a checkpoint file
    if (level == 0 && write_ascii == 1) {
      // TODO: Would be nice to be able to use file_name_digits
      // instead of doing this
      int strlen = dir.length();
      // Remove the ".temp" from the directory
      std::string dirout = dir.substr(0, strlen - 5);
      size_t num_start_loc = dirout.find_last_not_of("0123456789") + 1;
      std::string dir_path = dir;
      size_t num_end_path = dir_path.find_last_of("/") + 1;
      dir_path = dir_path.substr(0, num_end_path);
      std::string fname =
        dir_path + "spray" + dirout.substr(num_start_loc, strlen) + ".p3d";
      WriteAsciiFile(fname);
    }
  }

  ///
  /// Set the value of particle state for all partilces on a level
  ///
  // void SetAll (amrex::Real val, int pstate_idx, int lev);

  ///
  /// Compute a maximum time step based on the particle velocities and a
  /// particle CFL number
  ///
  amrex::Real estTimestep(int level, amrex::Real cfl) const;

  ///
  /// Set appropriate number of ghost cells
  ///
  inline void setSprayGridInfo(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const int amr_iteration,
    const int fine_ref_ratio,
    int& ghost_width,
    int& where_width,
    int& state_width,
    int& source_width)
  {
    int stencil_deposition_width = 1;
    int stencil_interpolation_width = 1;
    // Set the number of cells on level+1 to hold ghost particles
    ghost_width = 0;
    if (level < finest_level)
      ghost_width = fine_ref_ratio;
    int ghost_num = 0;
    if (level > 0) {
      ghost_num += amr_ncycle + amr_ncycle / 2;
      // For the initial iters, this will not work so we much have a check
      ghost_num = amrex::max(3, ghost_num);
    }
    // Set the number of ghost cells needed for interpolation
    state_width = stencil_interpolation_width + ghost_num;
    // Set the number of ghost cells needed for source deposition
    source_width = stencil_deposition_width + ghost_num;
    // Set the number of cells used in the where call
    // TODO: This is now unused. Only caused issues
    where_width = state_width;
  }

  ///
  /// Update the particles but do not move particles
  ///
  void moveKick(
    amrex::MultiFab& state,
    amrex::MultiFab& source,
    const int level,
    const amrex::Real& dt,
    const amrex::Real time,
    const bool isVirtualPart,
    const bool isGhostPart,
    const int state_ghosts,
    const int source_ghosts,
    amrex::MultiFab* u_mac = nullptr);

  ///
  /// Update and move positions of the particles
  ///
  void moveKickDrift(
    amrex::MultiFab& state,
    amrex::MultiFab& source,
    const int level,
    const amrex::Real& dt,
    const amrex::Real time,
    const bool isVirtualPart,
    const bool isGhostPart,
    const int state_ghosts,
    const int source_ghosts,
    const bool do_move,
    const int where_width,
    amrex::MultiFab* u_mac = nullptr);

  ///
  /// Update particles
  ///
  void updateParticles(
    const int& level,
    amrex::MultiFab& state,
    amrex::MultiFab& source,
    const amrex::Real& flow_dt,
    const amrex::Real& time,
    const int state_ghosts,
    const int source_ghosts,
    const bool isActive,
    const bool do_move,
    amrex::MultiFab* u_mac);

  // Modify particles based on walls
  // This creates new particles from splashing,
  // consolidates particles to wall films,
  // and removes particles that have left the domain
  void wallImpingement(
    const int& level,
    const amrex::Real& flow_dt,
    const amrex::Real& time,
#ifdef AMREX_USE_EB
    const amrex::FabArray<amrex::EBCellFlagFab>& flagmf,
    const amrex::MultiCutFab* bndrycent,
    const amrex::MultiCutFab* bndrynorm,
#endif
    const int state_ghosts,
    const int source_ghosts,
    const bool isActive);

  ///
  /// Sum boundary and copy temporary source term
  ///
  inline void transferSource(
    const int& source_ghosts,
    const int& level,
    amrex::MultiFab& tmp_source,
    amrex::MultiFab& source)
  {
    BL_PROFILE("ParticleContainer::transferSource()");
    const int nghost = amrex::min(source.nGrow(), source_ghosts);
    const int ncomp = tmp_source.nComp();
#ifdef SPRAY_PELE_LM
    if (level > 0) {
      amrex::MultiFab tmp_mf(
        this->m_gdb->ParticleBoxArray(level),
        this->m_gdb->ParticleDistributionMap(level), ncomp, source_ghosts);
      tmp_mf.setVal(0.);
      tmp_mf.copy(
        tmp_source, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity(), FabArrayBase::ADD);
      tmp_source.copy(
        tmp_mf, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity(), FabArrayBase::COPY);
    } else
#endif
    {
      tmp_source.SumBoundary(Geom(level).periodicity());
    }
    amrex::MultiFab::Add(source, tmp_source, 0, 0, source.nComp(), nghost);
    tmp_source.setVal(0.);
  }

  bool injectParticles(
    amrex::Real time,
    amrex::Real dt,
    int nstep,
    int level,
    int finest_level,
#ifdef SPRAY_PELE_LM
    ProbParm const& prob_parm
#else
    ProbParmHost const& prob_parm,
    ProbParmDevice const& prob_parm_d
#endif
  );

  void InitSprayParticles(
#ifdef SPRAY_PELE_LM
    ProbParm const& prob_parm
#else
    ProbParmHost const& prob_parm, ProbParmDevice const& prob_parm_d
#endif
  );

private:
  amrex::Real m_injectVel;
  // The number of spray droplets per computational particle
  amrex::Real m_parcelSize;
  // Temperature of walls
  amrex::Real m_wallT;
  ///
  /// This defines reflect_lo and reflect_hi from phys_bc
  ///
  void init_bcs();

  amrex::BCRec* phys_bc;
  bool reflect_lo[AMREX_SPACEDIM];
  bool reflect_hi[AMREX_SPACEDIM];
  SprayData* m_sprayData;
  SprayData* d_sprayData;
  SprayComps m_sprayIndx;
};

#endif
