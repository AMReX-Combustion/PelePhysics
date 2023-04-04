#ifndef SPRAYPARTICLES_H
#define SPRAYPARTICLES_H

#include "SprayFuelData.H"
#include <AMReX_AmrParticles.H>
#include <AMReX_Geometry.H>
#include "SprayJet.H"

#ifdef PELELM_USE_SPRAY
#include "pelelm_prob_parm.H"
#else
#include "prob_parm.H"
#endif

// Need components for velocity, diameter, temperature, mass fractions,
// breakup model variables, and wall film volume
#define NSR_SPR (SprayComps::pstateNum)
#define NSI_SPR 0
#define NAR_SPR 0
#define NAI_SPR 0

// Forward declarations
class SBPtrs;

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
    amrex::Real part_cfl = -1.)
    : amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>(amr),
      phys_bc(_phys_bc),
      m_sprayData(new SprayData{}),
      m_sprayIndx(SPI),
      m_partCFL(part_cfl)
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

  /// \brief Generalized injection routine for a single SprayJet
  void sprayInjection(
    const amrex::Real time,
    SprayJet* spray_jet,
    const amrex::Real sim_dt,
    const int level);

  /// \brief General initialization routine for uniformly distributed droplets
  void uniformSprayInit(
    const amrex::IntVect num_part,
    const amrex::RealVect vel_part,
    const amrex::Real dia_part,
    const amrex::Real T_part,
    const amrex::Real* Y_part,
    const int level,
    const int numRedist,
    const amrex::Real num_ppp = 1.);

  /// \brief Setup spray parameters
  static void spraySetup(SprayData& sprayData, const amrex::Real* body_force);

  /// \brief Read in spray parameters from input file
  static void readSprayParams(
    int& particle_verbose,
    amrex::Real& particle_cfl,
    int& write_spray_ascii_files,
    int& plot_spray_src,
    int& init_function,
    std::string& init_file,
    SprayData& sprayData,
    const amrex::Real& max_cfl = 0.5);

  /// \brief Spray particle write routine, writes plot, checkpoint, ascii, and
  /// injection data files
  void SprayParticleIO(
    const int level,
    const bool is_checkpoint,
    const int write_ascii,
    const std::string& dir);

  /// \brief Derive grid variables related to sprays
  void computeDerivedVars(
    amrex::MultiFab& mf_var, const int level, const int start_indx);

  /// \brief Compute a maximum time step based on the particle velocities and a
  /// particle CFL number
  amrex::Real estTimestep(int level, amrex::Real cfl) const;

  /// \brief Reset the particle ID in case we need to reinitialize the particles
  static inline void resetID(const int id) { ParticleType::NextID(id); }

  /// \brief Returns the number of ghost cells for making ghost particles. This
  /// is called on level N-1 to make ghost particles on level N from valid
  /// particles on level N-1
  /// @param level Level N
  /// @param finest_level Maximum level
  /// @param amr_ncycle Number of AMR subcycles taken on the Eulerian mesh
  /// @param cfl Particle CFL number on level N
  /// @param depos_width Number of cells adjacent to the cell containing the
  /// particle needed to interpolate the particle source term to the mesh
  static inline int getGhostPartCells(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const amrex::Real& cfl = 0.5,
    const int depos_width = 0)
  {
    if (level - 1 < finest_level) {
      return depos_width + static_cast<int>(std::round(
                             cfl * static_cast<amrex::Real>(amr_ncycle)));
    }
    return 0;
  }

  /// \brief Returns the number of ghost cells needed for a MultiFab used to
  /// interpolate a state to a particle location
  /// @param level Level N
  /// @param finest_level Maximum level
  /// @param amr_ncycle Number of AMR subcycles taken on the Eulerian mesh
  /// @param cfl Particle CFL number on level N
  /// @param interp_width Number of cells adjacent to the cell containing the
  /// particle needed to interpolated the state to the particle location; this
  /// is 1 for trilinear interpolation
  /// @param depos_width Number of cells adjacent to the cell containing the
  /// particle needed to interpolate the particle source term to the mesh
  static inline int getStateGhostCells(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const amrex::Real& cfl = 0.5,
    const int interp_width = 1,
    const int depos_width = 0)
  {
    int ghost_state = interp_width + static_cast<int>(std::round(cfl));
    if (level > 0) {
      // If ghost particles are present, we need to accommodate those
      const int ghost_part_cells =
        getGhostPartCells(level, finest_level, amr_ncycle, cfl, depos_width);
      int ghost_part_state = interp_width + ghost_part_cells;
      ghost_state = amrex::max(ghost_part_state, ghost_state);
    }
    return ghost_state;
  }

  /// \brief Returns the number of ghost cells needed for a particle to
  /// distribute the source terms to a Eulerian mesh
  /// @param level Level N
  /// @param finest_level Maximum level
  /// @param amr_ncycle Number of AMR subcycles taken on the Eulerian mesh
  /// @param cfl Particle CFL number on level N
  /// @param depos_width Number of cells adjacent to the cell containing the
  /// particle needed to interpolate the particle source term to the mesh
  static inline int getSourceGhostCells(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const amrex::Real& cfl = 0.5,
    const int depos_width = 0)
  {
    int ghost_source = depos_width + static_cast<int>(std::round(cfl));
    if (level > 0) {
      // If ghost particles are present, we need to accommodate those
      const int ghost_part_cells =
        getGhostPartCells(level, finest_level, amr_ncycle, cfl, depos_width);
      ghost_source = amrex::max(ghost_part_cells, ghost_source);
    }
    return ghost_source;
  }

  /// \brief Update the particles but do not move particles
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

  /// \brief Update and move positions of the particles
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

  /// \brief Update particle data and solve gas-phase source terms
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

  /// \brief Sum boundary and copy temporary source term
  inline void transferSource(
    const int& source_ghosts,
    const int& level,
    amrex::MultiFab& sprayDot,
    amrex::MultiFab& extSrc)
  {
    BL_PROFILE("ParticleContainer::transferSource()");
    AMREX_ALWAYS_ASSERT(level >= 0);
    const int nghost = amrex::min(extSrc.nGrow(), source_ghosts);
#ifdef PELELM_USE_SPRAY
    const int ncomp = sprayDot.nComp();
    if (level > 0) {
      amrex::MultiFab tmp_mf(
        this->m_gdb->ParticleBoxArray(level),
        this->m_gdb->ParticleDistributionMap(level), ncomp, source_ghosts);
      tmp_mf.setVal(0.);
      tmp_mf.ParallelAdd(
        sprayDot, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity());
      sprayDot.ParallelCopy(
        tmp_mf, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity());
    } else
#endif
    {
      sprayDot.SumBoundary(Geom(level).periodicity());
    }
    if (sprayDot.nComp() == extSrc.nComp()) {
      amrex::MultiFab::Add(
        extSrc, sprayDot, 0, 0, extSrc.nComp(), nghost); // NOLINT
    } else {
      for (amrex::MFIter mfi(extSrc, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        const amrex::Box vbox = mfi.growntilebox(nghost);
        auto const& spraydot = sprayDot.array(mfi);
        auto const& ext_src = extSrc.array(mfi);
        addSpraySrc(vbox, spraydot, ext_src);
      }
    }
    sprayDot.setVal(0.);
  }

  /// \brief Add spray source term data with components pstateNum to data
  /// containing total number of conservative variables
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
          const int ecomp = SPI.specIndx + fdat->dep_indx[sp];
          ext_src(i, j, k, ecomp) += spraydot(i, j, k, SPI.specSrcIndx + sp);
        }
      });
  }

  /// \brief Problem specific injection routine, calls #sprayInjection
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

  /// \brief Problem specific initialization routine. Typically initializes the
  /// SprayJet pointers or calls #uniformSprayInit
  void InitSprayParticles(
    const bool init_parts,
#ifdef PELELM_USE_SPRAY
    ProbParm const& prob_parm
#else
    ProbParmHost const& prob_parm,
    ProbParmDevice const& prob_parm_d
#endif
  );

  /// \brief Should be called after Restart or initialize routine. Reads
  /// injection data files if they are present. Checks to ensure all jet names
  /// are unique
  /// @param dir Name of restart directory
  void PostInitRestart(const std::string& dir = "");

  static std::string spray_fuel_names[SPRAY_FUEL_NUM];
  static std::string spray_dep_names[SPRAY_FUEL_NUM];
  static amrex::Vector<std::string> spray_derive_vars;
  static amrex::Real max_num_ppp;
  // The number density of a parcel formed from breakup is equal to
  // breakup_ppp_fact * parent parcel number density
  static amrex::Real breakup_ppp_fact;
  static amrex::Real B0_KHRT;
  static amrex::Real B1_KHRT;
  static amrex::Real C3_KHRT;

private:
  /// \brief This defines reflect_lo and reflect_hi from phys_bc
  void init_bcs();

  amrex::BCRec* phys_bc;
  bool reflect_lo[AMREX_SPACEDIM];
  bool reflect_hi[AMREX_SPACEDIM];
  SprayData* m_sprayData;
  SprayData* d_sprayData = nullptr;
  SprayComps m_sprayIndx;
  // Particle CFL, used for injection cases with LMeX
  amrex::Real m_partCFL;
  amrex::Real m_injectVel = 0.;
  amrex::Vector<std::unique_ptr<SprayJet>> m_sprayJets;
};

#endif
