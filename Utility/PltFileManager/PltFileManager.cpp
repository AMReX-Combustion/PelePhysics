#include <PltFileManager.H>
#include <AMReX_VisMF.H>
#include <PltFileManagerBCFill.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <utility>

using namespace amrex;

namespace pele {
namespace physics {
namespace pltfilemanager {

namespace {
const std::string level_prefix{"Level_"};
}

void
GotoNextLine(std::istream& is)
{
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}

PltFileManager::PltFileManager(std::string a_pltFile)
  : m_pltFile{std::move(a_pltFile)}
{
  // Get the plt metadata and resize part of the data vectors
  std::string pltFileHeader(m_pltFile + "/Header");
  readGenericPlotfileHeader(pltFileHeader);

  // Resize the actual data container
  m_dmaps.resize(m_nlevels);
  m_data.resize(m_nlevels);

  // Read the pltfile data
  readPlotFile();
}

void
PltFileManager::readGenericPlotfileHeader(const std::string& a_pltFileHeader)
{
  Vector<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(a_pltFileHeader, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);

  std::string line, word;

  // Title line
  std::getline(is, line);

  // Number of variables
  m_nvars = 0;
  is >> m_nvars;
  GotoNextLine(is);

  // Extract variables names
  m_vars.resize(m_nvars);
  for (int n = 0; n < m_nvars; n++) {
    is >> m_vars[n];
    GotoNextLine(is);
  }

  // Get and check space dimension
  int PLT_SPACEDIM = AMREX_SPACEDIM;
  is >> PLT_SPACEDIM;
  GotoNextLine(is);
  AMREX_ASSERT(PLT_SPACEDIM == AMREX_SPACEDIM);

  // Simulation time
  is >> m_time;
  GotoNextLine(is);

  // Number of levels
  is >> m_nlevels;
  GotoNextLine(is);
  m_nlevels += 1; // Finest is stored, need to add 1

  // Setup data holders
  m_grids.resize(m_nlevels);
  m_geoms.resize(m_nlevels);
  m_refRatio.resize(m_nlevels - 1);

  // Level 0 geometry
  Real prob_lo[AMREX_SPACEDIM];
  Real prob_hi[AMREX_SPACEDIM];
  // Low coordinates of domain bounding box
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      prob_lo[i++] = std::stod(word);
    }
  }

  // High coordinates of domain bounding box
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      prob_hi[i++] = std::stod(word);
    }
  }

  // Set up PltFile domain real box
  RealBox rb(prob_lo, prob_hi);

  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      m_refRatio[i++] = std::stoi(word);
    }
  }

  // Get levels Domains
  Vector<Box> Domains(m_nlevels);
  for (int lev = 0; lev < m_nlevels; ++lev) {
    is >> Domains[lev];
  }
  GotoNextLine(is);
  GotoNextLine(is); // Skip nsteps line
  for (int lev = 0; lev < m_nlevels; ++lev) {
    GotoNextLine(is); // Skip dx line
  }

  // Coordinate system
  int coord_sys = 0;
  is >> coord_sys;
  GotoNextLine(is);

  // Populate the geometry vector, assume no periodicity
  Array<int, AMREX_SPACEDIM> perio({AMREX_D_DECL(0, 0, 0)});
  m_geoms[0] = Geometry(Domains[0], rb, coord_sys, perio);
  for (int lev = 1; lev < m_nlevels; ++lev) {
    m_geoms[lev] = refine(m_geoms[lev - 1], m_refRatio[lev - 1]);
  }
}

void
PltFileManager::readPlotFile()
{
  // Set BoxArray, DistMap on each level and load data
  // TODO: only load a subset of the pltfile variables
  for (int lev = 0; lev < m_nlevels; ++lev) {
    readLevelBoxArray(lev, m_grids[lev]);
    m_dmaps[lev] = DistributionMapping(m_grids[lev]);
    m_data[lev].define(m_grids[lev], m_dmaps[lev], m_nvars, 0);
    VisMF::Read(
      m_data[lev],
      MultiFabFileFullPrefix(lev, m_pltFile, level_prefix, "Cell"));
  }
}

void
PltFileManager::readLevelBoxArray(int a_lev, BoxArray& a_grid)
{
  const std::string lvlHeader(
    m_pltFile + "/" + level_prefix + std::to_string(a_lev) + "/Cell_H");

  Vector<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(lvlHeader, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);

  std::string line;

  std::getline(is, line); // Skip
  GotoNextLine(is);       // Skip
  GotoNextLine(is);       // Skip
  GotoNextLine(is);       // Skip
  a_grid.readFrom(is);
}

void
PltFileManager::fillPatchFromPlt(
  int a_lev,
  const Geometry& a_level_geom,
  int pltComp,
  int dataComp,
  int nComp,
  MultiFab& a_mf)
{
  Vector<BCRec> dummyBCRec(nComp);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    if (a_level_geom.isPeriodic(idim)) {
      for (int n = 0; n < nComp; n++) {
        dummyBCRec[n].setLo(idim, BCType::int_dir);
        dummyBCRec[n].setHi(idim, BCType::int_dir);
      }
    } else {
      for (int n = 0; n < nComp; n++) {
        dummyBCRec[n].setLo(idim, BCType::foextrap);
        dummyBCRec[n].setHi(idim, BCType::foextrap);
      }
    }
  }

  // There might be a number of problems related to proper nesting of the
  // current BA and the PltFile BA. Try to address some of those, but need
  // further testing.

  if (a_lev == 0) {
    // Ensure that the provided level 0 geometry is contained in the one
    // we've read from the pltfile
    AMREX_ALWAYS_ASSERT(
      m_geoms[0].ProbDomain().contains(a_level_geom.ProbDomain(), 0.0000001));

    // Check the refRatio between PltFile level 0 and ours
    IntVect lev0rr = a_level_geom.Domain().size() / m_geoms[0].Domain().size();

    // Same domain size, just do a fill patch single level
    if (lev0rr == IntVect::TheUnitVector()) {
      PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> bndry_func(
        a_level_geom, {dummyBCRec}, FillExtDirDummy{});
      FillPatchSingleLevel(
        a_mf, IntVect(0), 0.0, {&(m_data[a_lev])}, {0.0}, pltComp, dataComp,
        nComp, a_level_geom, bndry_func, 0);

      // Our level 0 is finer than the PltFile one.
    } else if (lev0rr.max() > 1) {

      // Interpolator (need EB version ?)
      InterpBase* mapper = &mf_cell_cons_interp;

      // Start by filling all the data with a coarseInterp. We've checked that
      // the geom RealBoxes match already, so PltFile level 0 is good to interp.
      {
        PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> crse_bndry_func(
          m_geoms[0], {dummyBCRec}, FillExtDirDummy{});
        PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> fine_bndry_func(
          a_level_geom, {dummyBCRec}, FillExtDirDummy{});
        InterpFromCoarseLevel(
          a_mf, IntVect(0), 0.0, m_data[0], pltComp, dataComp, nComp,
          m_geoms[0], a_level_geom, crse_bndry_func, 0, fine_bndry_func, 0,
          lev0rr, mapper, {dummyBCRec}, 0);
      }

      // Then get data from the PltFile finer levels if any
      for (int pltlev = 1; pltlev < m_geoms.size(); pltlev++) {

        IntVect rr =
          a_level_geom.Domain().size() / m_geoms[pltlev].Domain().size();

        // Current Plt level resolution matches our: lets interp and wrap it up
        if (rr == IntVect::TheUnitVector()) {
          MultiFab temp(m_grids[pltlev], m_dmaps[pltlev], nComp, 0);
          PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> bndry_func(
            m_geoms[pltlev], {dummyBCRec}, FillExtDirDummy{});
          FillPatchSingleLevel(
            temp, IntVect(0), 0.0, {&m_data[pltlev]}, {0.0}, pltComp, 0, nComp,
            m_geoms[pltlev], bndry_func, 0);

          a_mf.ParallelCopy(temp, 0, dataComp, nComp);

          break;
        }

        // Otherwise do another InterpFromCoarseLevel
        MultiFab temp(refine(m_grids[pltlev], rr), m_dmaps[pltlev], nComp, 0);
        PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> crse_bndry_func(
          m_geoms[pltlev], {dummyBCRec}, FillExtDirDummy{});
        PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> fine_bndry_func(
          a_level_geom, {dummyBCRec}, FillExtDirDummy{});
        InterpFromCoarseLevel(
          temp, IntVect(0), 0.0, m_data[pltlev], pltComp, 0, nComp,
          m_geoms[pltlev], a_level_geom, crse_bndry_func, 0, fine_bndry_func, 0,
          rr, mapper, {dummyBCRec}, 0);
        a_mf.ParallelCopy(temp, 0, dataComp, nComp);
      }

    } else {
      // Otherwise bail out because we don't handle averaging_down (TODO ?)
      Abort("When initializing dataFromPlt, our level 0 can't be coarser than "
            "the PltFile one");
    }
  } else {
    // Check the refRatio between PltFile level 0 and the current level
    IntVect lev0rr = a_level_geom.Domain().size() / m_geoms[0].Domain().size();

    // Interpolator (need EB version ?)
    InterpBase* mapper = &mf_cell_cons_interp;

    // Start by filling the entire level with level 0 from PltFile, to ensure we
    // have some data everywhere
    {
      PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> crse_bndry_func(
        m_geoms[0], {dummyBCRec}, FillExtDirDummy{});
      PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> fine_bndry_func(
        a_level_geom, {dummyBCRec}, FillExtDirDummy{});
      InterpFromCoarseLevel(
        a_mf, IntVect(0), 0.0, m_data[0], pltComp, dataComp, nComp, m_geoms[0],
        a_level_geom, crse_bndry_func, 0, fine_bndry_func, 0, lev0rr, mapper,
        {dummyBCRec}, 0);
    }

    // Then get data from the PltFile finer levels if any
    for (int pltlev = 1; pltlev < m_geoms.size(); pltlev++) {

      IntVect rr =
        a_level_geom.Domain().size() / m_geoms[pltlev].Domain().size();

      // Current Plt level resolution matches our: lets interp and wrap it up
      if (rr == IntVect::TheUnitVector()) {
        MultiFab temp(m_grids[pltlev], m_dmaps[pltlev], nComp, 0);
        PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> bndry_func(
          m_geoms[pltlev], {dummyBCRec}, FillExtDirDummy{});
        FillPatchSingleLevel(
          temp, IntVect(0), 0.0, {&m_data[pltlev]}, {0.0}, pltComp, 0, nComp,
          m_geoms[pltlev], bndry_func, 0);

        a_mf.ParallelCopy(temp, 0, dataComp, nComp);

        break;
      }

      // Otherwise do another InterpFromCoarseLevel
      MultiFab temp(refine(m_grids[pltlev], rr), m_dmaps[pltlev], nComp, 0);
      PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> crse_bndry_func(
        m_geoms[pltlev], {dummyBCRec}, FillExtDirDummy{});
      PhysBCFunct<GpuBndryFuncFab<FillExtDirDummy>> fine_bndry_func(
        a_level_geom, {dummyBCRec}, FillExtDirDummy{});
      InterpFromCoarseLevel(
        temp, IntVect(0), 0.0, m_data[pltlev], pltComp, 0, nComp,
        m_geoms[pltlev], a_level_geom, crse_bndry_func, 0, fine_bndry_func, 0,
        rr, mapper, {dummyBCRec}, 0);
      a_mf.ParallelCopy(temp, 0, dataComp, nComp);
    }
  }
}

} // namespace pltfilemanager
} // namespace physics
} // namespace pele
