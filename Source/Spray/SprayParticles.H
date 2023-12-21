#ifndef SPRAYPARTICLES_H
#define SPRAYPARTICLES_H

#include "SprayFuelData.H"
#include <AMReX_AmrParticles.H>
#include <AMReX_Geometry.H>
#include "SprayJet.H"

// Need components for velocity, diameter, temperature, mass fractions,
// breakup model variables, and wall film volume
#define NSR_SPR (SprayComps::pstateNum)
#define NSI_SPR 0
#define NAR_SPR 0
#define NAI_SPR 0

// Forward declarations
struct SBPtrs;

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

  SprayParticleContainer(amrex::AmrCore* amr, amrex::BCRec* _phys_bc)
    : amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>(amr),
      phys_bc(_phys_bc)
  {
    init_bcs();
  }

  static void SprayCleanUp()
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
  /// @param num_part Number of parcels to initialize in each direction
  /// @param vel_part Droplet velocity
  /// @param dia_part Droplet diameter
  /// @param T_part Droplet temperature
  /// @param Y_part Pointer to array of droplet mass fractions
  /// @param level Current AMR level
  /// @param num_redist Number of redistributes to do during initialization
  /// @param num_ppp Number of droplets per parcel
  void uniformSprayInit(
    const amrex::IntVect num_part,
    const amrex::RealVect vel_part,
    const amrex::Real dia_part,
    const amrex::Real T_part,
    const amrex::Real* Y_part,
    const int level,
    const int num_redist = 1,
    const amrex::Real num_ppp = 1.);

  /// \brief Setup spray parameters
  static void spraySetup(const amrex::Real* body_force);

  /// \brief Read in spray parameters from input file
  static void readSprayParams(int& particle_verbose);

  /// \brief Create droplets from splashing or breakup
  void CreateSBDroplets(
    const int Np,
    const amrex::Real sub_dt,
    const splash_breakup* N_SB_h,
    const SBPtrs& rfh,
    const int level);

  /// \brief Spray particle write routine, writes plot, checkpoint, ascii, and
  /// injection data files
  void SprayParticleIO(
    const int level, const bool is_checkpoint, const std::string& dir);

  /// \brief Derive grid variables related to sprays
  void computeDerivedVars(
    amrex::MultiFab& mf_var, const int level, const int start_indx);

  /// \brief Compute a maximum time step based on the particle velocities and a
  /// particle CFL number
  amrex::Real estTimestep(int level) const;

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
      return amrex::max(
        1, depos_width + static_cast<int>(std::round(
                           cfl * static_cast<amrex::Real>(amr_ncycle))));
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
    int ghost_source =
      amrex::max(1, depos_width + static_cast<int>(std::round(cfl)));
    if (level > 0) {
      // If ghost particles are present, we need to accommodate those
      const int ghost_part_cells =
        getGhostPartCells(level, finest_level, amr_ncycle, cfl, depos_width);
      ghost_source = amrex::max(ghost_part_cells, ghost_source);
    }
    return ghost_source;
  }

  /// \brief Update but do not move particles
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
    const amrex::Real spray_cfl_lev = -1.);

  /// \brief Update and move particles
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
    const amrex::Real spray_cfl_lev = -1.);

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
    const amrex::Real spray_cfl_lev);

  /// \brief Sum boundary and copy temporary source term
  /// @param source_ghosts Number of ghost cells needed for source terms
  /// @param level Current AMR level
  /// @param tmpSource MultiFab containing the temporary spray source data
  /// @param actSource MultiFab where final source data should end up
  inline void transferSource(
    const int& source_ghosts,
    const int& level,
    amrex::MultiFab& tmpSource,
    amrex::MultiFab& actSource)
  {
    BL_PROFILE("ParticleContainer::transferSource()");
    AMREX_ALWAYS_ASSERT(level >= 0);
    const int nghost = amrex::min(actSource.nGrow(), source_ghosts);
#ifdef PELELM_USE_SPRAY
    const int ncomp = tmpSource.nComp();
    if (level > 0) {
      amrex::MultiFab tmp_mf(
        this->m_gdb->ParticleBoxArray(level),
        this->m_gdb->ParticleDistributionMap(level), ncomp, source_ghosts);
      tmp_mf.setVal(0.);
      tmp_mf.ParallelAdd(
        tmpSource, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity());
      tmpSource.ParallelCopy(
        tmp_mf, 0, 0, ncomp, source_ghosts, source_ghosts,
        Geom(level).periodicity());
    } else
#endif
    {
      tmpSource.SumBoundary(Geom(level).periodicity());
    }
    if (tmpSource.nComp() == actSource.nComp()) {
      amrex::MultiFab::Add(
        actSource, tmpSource, 0, 0, actSource.nComp(), nghost); // NOLINT
    } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        for (amrex::MFIter mfi(actSource, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
          const amrex::Box vbox = mfi.growntilebox(nghost);
          auto const& tmp_array = tmpSource.array(mfi);
          auto const& act_array = actSource.array(mfi);
          addSpraySrc(vbox, tmp_array, act_array);
        }
      }
    }
    tmpSource.setVal(0.);
  }

  /// \brief Add spray source term data with components pstateNum to data
  /// containing total number of conservative variables
  static void addSpraySrc(
    const amrex::Box& vbox,
    amrex::Array4<const amrex::Real> const& tmp_array,
    amrex::Array4<amrex::Real> const& act_array)
  {
    SprayData* fdat = d_sprayData;
    SprayComps SPI = m_sprayIndx;
    amrex::ParallelFor(
      vbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          act_array(i, j, k, SPI.momIndx + dir) +=
            tmp_array(i, j, k, SPI.momSrcIndx + dir);
        }
        act_array(i, j, k, SPI.rhoIndx) += tmp_array(i, j, k, SPI.rhoSrcIndx);
        act_array(i, j, k, SPI.engIndx) += tmp_array(i, j, k, SPI.engSrcIndx);
        for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
          const int ecomp = SPI.specIndx + fdat->dep_indx[sp];
          act_array(i, j, k, ecomp) += tmp_array(i, j, k, SPI.specSrcIndx + sp);
        }
      });
  }

  /// \brief Problem specific injection routine, calls #sprayInjection
  static bool injectParticles(
    amrex::Real time, amrex::Real dt, int nstep, int level, int finest_level);

  /// \brief Generalized initializer routine. Calls #InitSprayParticles and
  /// #PostInitRestart
  void SprayInitialize(const std::string& restart_dir = "");

  /// \brief Problem specific initialization routine. Typically initializes the
  /// SprayJet pointers or calls #uniformSprayInit
  void InitSprayParticles(const bool init_parts);

  /// \brief Should be called after Restart or initialize routine. Reads
  /// injection data files if they are present. Checks to ensure all jet names
  /// are unique
  /// @param dir Name of restart directory
  void PostInitRestart(const std::string& dir = "");

  /// \brief Return the number of derived spray variables
  static int NumDeriveVars()
  {
    return static_cast<int>(m_sprayDeriveVars.size());
  }

  /// \brief Return derive variable names
  static amrex::Vector<std::string>& DeriveVarNames()
  {
    return m_sprayDeriveVars;
  }

  /// \brief Return host spray data pointer
  static SprayData* getSprayData() { return m_sprayData; }

  /// \brief Return fuel species index
  static int getFuelIndx(const int spf) { return m_sprayData->dep_indx[spf]; }

  /// \brief Return spray indices
  static SprayComps getSprayComps() { return m_sprayIndx; }

  static void AssignSprayComps(SprayComps scomps) { m_sprayIndx = scomps; }

  static std::string m_sprayFuelNames[SPRAY_FUEL_NUM];
  static std::string m_sprayDepNames[SPRAY_FUEL_NUM];
  static amrex::Vector<std::string> m_sprayDeriveVars;
  static amrex::Real m_maxNumPPP;
  // Factor for balancing number of child parcels formed versus child parcel
  // number density 1: Fewer parcels, higher number densities 0: More parcels,
  // lower number densities
  static amrex::Real m_breakupPPPFact;
  static amrex::Real m_khrtB0;
  static amrex::Real m_khrtB1;
  static amrex::Real m_khrtC3;
  static SprayData* m_sprayData;
  static SprayData* d_sprayData;
  static SprayComps m_sprayIndx;
  static amrex::Real spray_cfl;
  static bool write_ascii_files;
  static bool plot_spray_src;
  static std::string spray_init_file;

private:
  /// \brief This defines reflect_lo and reflect_hi from phys_bc
  void init_bcs();

  amrex::BCRec* phys_bc;
  bool reflect_lo[AMREX_SPACEDIM];
  bool reflect_hi[AMREX_SPACEDIM];
  amrex::Vector<std::unique_ptr<SprayJet>> m_sprayJets;
};

#endif
