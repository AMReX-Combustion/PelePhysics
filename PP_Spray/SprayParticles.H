#ifndef SPRAYPARTICLES_H
#define SPRAYPARTICLES_H

#include "SprayFuelData.H"
#include <AMReX_AmrParticles.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_Particles.H>
#include <memory>

#ifdef PELELM_USE_SPRAY
#include "pelelm_prob.H"
#include "pelelm_prob_parm.H"
#else
#include "prob.H"
#include "prob_parm.H"
#endif

#define NSR_SPR (AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM)
#define NSI_SPR 0
#define NAR_SPR 0
#define NAI_SPR 0

class MyParIter : public amrex::ParIter<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>
{
public:
  using amrex::ParIter<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>::ParIter;
};

class MyParConstIter
  : public amrex::ParConstIter<NSR_SPR, NSI_SPR, NAR_SPR, NSI_SPR>
{
public:
  using amrex::ParConstIter<NSR_SPR, NSI_SPR, NAR_SPR, NSI_SPR>::ParConstIter;
};

class SprayParticleContainer
  : public amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>
{
public:
  using PairIndex = std::pair<int, int>;
  using HostVectReal = amrex::Gpu::HostVector<amrex::Real>;
  using HostVectInt = amrex::Gpu::HostVector<int>;

  SprayParticleContainer(
    amrex::AmrCore* amr,
    amrex::BCRec* _phys_bc,
    const SprayData& fdat,
    const SprayComps& SPI,
    amrex::Real wall_temp,
    amrex::Real part_cfl = -1.)
    : amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>(amr),
      phys_bc(_phys_bc),
      m_sprayData(new SprayData{}),
      d_sprayData(nullptr),
      m_sprayIndx(SPI),
      m_wallT(wall_temp),
      m_partCFL(part_cfl),
      m_injectVel(0.)
  {
    d_sprayData =
      static_cast<SprayData*>(amrex::The_Arena()->alloc(sizeof(SprayData)));
    m_sprayData->build(fdat);
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, m_sprayData, m_sprayData + 1, d_sprayData);
    init_bcs();
  }

  ~SprayParticleContainer() override
  {
    delete m_sprayData;
    amrex::The_Arena()->free(d_sprayData);
  }

  template <class Dist>
  void sprayInjection(
    Dist& dist_type,
    const amrex::RealVect jet_cent,
    const amrex::RealVect jet_norm,
    const amrex::Real jet_dia,
    const amrex::Real part_temp,
    const amrex::Real mass_flow,
    const amrex::Real jet_vel,
    const amrex::Real spread_angle,
    const amrex::Real dt,
    const amrex::Real* Y_jet,
    const int level,
    int startInjProc = 0,
    int numInjProcs = 1,
    const bool hollow_spray = false);

  static void readSprayParams(
    int& particle_verbose,
    amrex::Real& particle_cfl,
    amrex::Real& wall_temp,
    bool& mass_trans,
    bool& mom_trans,
    int& write_spray_ascii_files,
    int& plot_spray_src,
    int& init_function,
    std::string& init_file,
    SprayData& sprayData,
    std::string* sprayFuelNames,
    amrex::Vector<std::string>& derivePlotVars,
    const amrex::Real& max_cfl =
      0.5); // Only overwrite this if particles are subcycled

  void SprayParticleIO(
    const int level,
    const bool is_checkpoint,
    const int write_ascii,
    const std::string& dir,
    const std::string* sprayFuelNames)
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
    for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
      real_comp_names[pstateY + sp] = "spray_mf_" + sprayFuelNames[sp];
    }
    amrex::Vector<std::string> int_comp_names;
    Checkpoint(
      dir, "particles", is_checkpoint, real_comp_names, int_comp_names);
    // Here we write ascii information every time we write a checkpoint file
    if (level == 0 && write_ascii == 1) {
      // TODO: Would be nice to be able to use file_name_digits
      // instead of doing this
      size_t num_end_loc = dir.find_last_of("0123456789") + 1;
      // Remove anything following numbers, like .temp
      std::string dirout = dir.substr(0, num_end_loc);
      size_t num_start_loc = dirout.find_last_not_of("0123456789") + 1;
      std::string numstring = dirout.substr(num_start_loc, dirout.length());
      std::string dir_path = dir;
      size_t num_end_path = dir_path.find_last_of("/") + 1;
      dir_path = dir_path.substr(0, num_end_path);
      std::string fname = dir_path + "spray" + numstring + ".p3d";
      WriteAsciiFile(fname);
    }
  }

  ///
  /// Derive grid variables related to sprays
  ///
  void computeDerivedVars(
    amrex::MultiFab& mf_var,
    const int level,
    const int start_indx,
    const amrex::Vector<std::string> derivePlotVars,
    const std::string* sprayFuelNames);

  ///
  /// Compute a maximum time step based on the particle velocities and a
  /// particle CFL number
  ///
  amrex::Real estTimestep(int level, amrex::Real cfl) const;

  ///
  /// Reset the particle ID in case we need to reinitialize the particles
  ///
  static inline void resetID(const int id) { ParticleType::NextID(id); }

  // This returns the number of ghost cells for making ghost particles
  // level - Level N
  // finest_level - Maximum level
  // amr_ncycle - Number of AMR subcycles taken on the Eulerian mesh
  // cfl - Particle CFL number on level N
  // depos_width - Number of source term cells influenced by a particle
  // NOTE: This is called on level N-1 to make ghost particles on level N
  // from valid particles on level N-1
  static inline int getGhostPartCells(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const amrex::Real& cfl = 0.5,
    const int depos_width = 1)
  {
    if (level - 1 < finest_level) {
      return depos_width + int(std::round(cfl * amrex::Real(amr_ncycle)));
    }
    return 0;
  }

  // This returns the number of ghost cells needed for a MultiFab
  // used to interpolate a state to a particle location
  static inline int getStateGhostCells(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const amrex::Real& cfl = 0.5,
    const int interp_width = 1,
    const int depos_width = 1)
  {
    int ghost_state = interp_width + int(std::round(cfl));
    if (level > 0) {
      // If ghost particles are present, we need to accommodate those
      const int ghost_part_cells =
        getGhostPartCells(level, finest_level, amr_ncycle, cfl, depos_width);
      int ghost_part_state = interp_width + ghost_part_cells;
      ghost_state = amrex::max(ghost_part_state, ghost_state);
    }
    return ghost_state;
  }

  // This returns the number of ghost cells needed for a MultiFab
  // used to interpolate a state to a particle location
  static inline int getSourceGhostCells(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const amrex::Real& cfl = 0.5,
    const int depos_width = 1)
  {
    int ghost_source = depos_width + int(std::round(cfl));
    if (level > 0) {
      // If ghost particles are present, we need to accommodate those
      const int ghost_part_cells =
        getGhostPartCells(level, finest_level, amr_ncycle, cfl, depos_width);
      ghost_source = amrex::max(ghost_part_cells, ghost_source);
    }
    return ghost_source;
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
    pele::physics::transport::TransParm<
      pele::physics::EosType,
      pele::physics::TransportType> const* ltransparm,
    const amrex::Real spray_cfl_lev = -1.,
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
    pele::physics::transport::TransParm<
      pele::physics::EosType,
      pele::physics::TransportType> const* ltransparm,
    const amrex::Real spray_cfl_lev = -1.,
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
    const bool isVirt,
    const bool isGhost,
    const bool do_move,
    pele::physics::transport::TransParm<
      pele::physics::EosType,
      pele::physics::TransportType> const* ltransparm,
    const amrex::Real spray_cfl_lev,
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
    AMREX_ALWAYS_ASSERT(level >= 0);
    const int nghost = amrex::min(source.nGrow(), source_ghosts);
#ifdef PELELM_USE_SPRAY
    const int ncomp = tmp_source.nComp();
    if (level > 0) {
      amrex::MultiFab tmp_mf(
        this->m_gdb->ParticleBoxArray(level),
        this->m_gdb->ParticleDistributionMap(level), ncomp, source_ghosts);
      tmp_mf.setVal(0.);
      tmp_mf.ParallelAdd(
        tmp_source, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity());
      tmp_source.ParallelCopy(
        tmp_mf, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity());
    } else
#endif
    {
      tmp_source.SumBoundary(Geom(level).periodicity());
    }
    if (tmp_source.nComp() == source.nComp()) {
      amrex::MultiFab::Add(source, tmp_source, 0, 0, source.nComp(), nghost);
    } else {
      for (amrex::MFIter mfi(source, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        const amrex::Box vbox = mfi.growntilebox(nghost);
        auto const& spraydot = tmp_source.array(mfi);
        auto const& ext_src = source.array(mfi);
        addSpraySrc(vbox, spraydot, ext_src);
      }
    }
    tmp_source.setVal(0.);
  }

  ///
  /// Add spray source term data (with components AMREX_SPACEDIM +
  /// SPRAY_FUEL_NUM + 2) to data containing total number of conservative
  /// variables
  ///
  inline void addSpraySrc(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& spraydot,
    amrex::Array4<amrex::Real> const& ext_src)
  {
    SprayData* fdat = d_sprayData;
    SprayComps SPI = m_sprayIndx;
    amrex::ParallelFor(
      vbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          ext_src(i, j, k, SPI.momIndx + dir) +=
            spraydot(i, j, k, SPI.momSrcIndx + dir);
        }
        ext_src(i, j, k, SPI.rhoIndx) += spraydot(i, j, k, SPI.rhoSrcIndx);
        ext_src(i, j, k, SPI.engIndx) += spraydot(i, j, k, SPI.engSrcIndx);
        for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
          const int ecomp = SPI.specIndx + fdat->indx[sp];
          ext_src(i, j, k, ecomp) += spraydot(i, j, k, SPI.specSrcIndx + sp);
        }
      });
  }

  bool injectParticles(
    amrex::Real time,
    amrex::Real dt,
    int nstep,
    int level,
    int finest_level,
#ifdef PELELM_USE_SPRAY
    ProbParm const& prob_parm
#else
    ProbParmHost const& prob_parm,
    ProbParmDevice const& prob_parm_d
#endif
  );

  void InitSprayParticles(
#ifdef PELELM_USE_SPRAY
    ProbParm const& prob_parm
#else
    ProbParmHost const& prob_parm, ProbParmDevice const& prob_parm_d
#endif
  );

private:
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
  // Temperature of walls
  amrex::Real m_wallT;
  // Particle CFL, used for injection cases with LMeX
  amrex::Real m_partCFL;
  amrex::Real m_injectVel;
};

#endif
