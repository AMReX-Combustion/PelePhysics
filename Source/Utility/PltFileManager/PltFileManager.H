#ifndef PLTFILEMANAGER_H
#define PLTFILEMANAGER_H

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>

namespace pele::physics::pltfilemanager {

class PltFileManager
{
public:
  explicit PltFileManager(std::string a_pltFile);

  void readGenericPlotfileHeader(const std::string& a_pltFileHeader);

  void readPlotFileMetaData();

  void readPlotFileData();

  void readLevelBoxArray(int a_lev, amrex::BoxArray& a_grid);

  void fillPatchFromPlt(
    int a_lev,
    const amrex::Geometry& a_level_geom,
    int pltComp,
    int dataComp,
    int ncomp,
    amrex::MultiFab& a_mf,
    int interp_type = 1);

  amrex::Vector<std::string> getVariableList() { return m_vars; }
  int getNlev() const { return m_nlevels; }
  amrex::Real getTime() const { return m_time; }
  int getNsteps() const { return m_nsteps; }
  const amrex::BoxArray& getGrid(int a_lev) { return m_grids[a_lev]; }
  const amrex::Geometry& getGeom(int a_lev) { return m_geoms[a_lev]; }
  int getRefRatio(int a_lev) { return m_refRatio[a_lev]; }

protected:
  // Plot metadata
  std::string m_pltFile;             // pltfile name
  int m_nvars{0};                    // number of variables
  amrex::Vector<std::string> m_vars; // list of variable in the plotfile
  int m_nlevels{0};                  // number of levels
  amrex::Real m_time{0.0};           // Simulation time
  int m_nsteps;                      // Number of steps
  bool m_dataLoaded{false}; // Flag to check is the data has been read in

  // Geometry, grid and data containers
  amrex::Vector<amrex::BoxArray> m_grids;
  amrex::Vector<amrex::Geometry> m_geoms;
  amrex::Vector<int> m_refRatio;
  amrex::Vector<amrex::DistributionMapping> m_dmaps;
  amrex::Vector<amrex::MultiFab> m_data;
};

} // namespace pele::physics::pltfilemanager
#endif
