#ifndef _SPRAYPARTICLES_H_
#define _SPRAYPARTICLES_H_

#include "SprayFuelData.H"
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

#define NSR_SPR AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM
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
    SprayData fdat,
    SprayComps SPI,
    amrex::Real wall_temp)
    : amrex::AmrParticleContainer<NSR_SPR, NSI_SPR, NAR_SPR, NAI_SPR>(amr),
      m_injectVel(0.),
      m_wallT(wall_temp),
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

  static void readSprayParams(
    int& particle_verbose,
    amrex::Real& particle_cfl,
    amrex::Real& wall_temp,
    int& mass_trans,
    int& mom_trans,
    int& write_spray_ascii_files,
    int& plot_spray_src,
    int& init_function,
    std::string& init_file,
    SprayData& sprayData,
    std::string* sprayFuelNames);

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
  /// Set the value of particle state for all partilces on a level
  ///
  // void SetAll (amrex::Real val, int pstate_idx, int lev);

  ///
  /// Compute a maximum time step based on the particle velocities and a
  /// particle CFL number
  ///
  amrex::Real estTimestep(int level, amrex::Real cfl) const;

  ///
  /// Reset the particle ID in case we need to reinitialize the particles
  ///
  inline void resetID(const int id) { ParticleType::NextID(id); }

  ///
  /// Set appropriate number of ghost cells
  ///
  static inline void setSprayGridInfo(
    const int level,
    const int finest_level,
    const int amr_ncycle,
    const int fine_ref_ratio,
    int& ghost_width,
    int& state_width,
    int& source_width)
  {
    int stencil_deposition_width = 1;
    int stencil_interpolation_width = 1;
    // Set the number of cells on level+1 to hold ghost particles
    ghost_width = 0;
    if (level < finest_level) {
      ghost_width = fine_ref_ratio;
    }
    int ghost_num = 0;
    if (level > 0) {
      ghost_num += amr_ncycle + amr_ncycle / 2;
    }
    // Set the number of ghost cells needed for interpolation
    state_width = stencil_interpolation_width + ghost_num;
    // Set the number of ghost cells needed for source deposition
    source_width = stencil_deposition_width + ghost_num;
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
    pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
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
#ifdef SPRAY_PELE_LM
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
