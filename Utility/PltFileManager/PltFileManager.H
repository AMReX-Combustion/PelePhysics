#ifndef _PLTFILEMANAGER_H_
#define _PLTFILEMANAGER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>

namespace pele {
namespace physics {
namespace pltfilemanager {

class PltFileManager
{
public:
  PltFileManager(const std::string& a_pltFile);

  void readGenericPlotfileHeader(const std::string& a_pltFileHeader);

  void readPlotFile();
  void readLevelBoxArray(int a_lev, amrex::BoxArray& a_grid);

  void fillPatchFromPlt(
    int a_lev,
    const amrex::Geometry& a_level_geom,
    int pltComp,
    int dataComp,
    int ncomp,
    amrex::MultiFab& a_mf);

  amrex::Vector<std::string> getVariableList() { return m_vars; }

protected:
  // MetaData
  std::string m_pltFile;             // pltfile name
  int m_nvars;                       // number of variables
  amrex::Vector<std::string> m_vars; // list of variable in the plotfile
  int m_nlevels;                     // number of levels
  amrex::Real m_time;                // Simulation time
  amrex::Real m_nsteps;              // Simulation step count

  amrex::Vector<amrex::BoxArray> m_grids;
  amrex::Vector<amrex::Geometry> m_geoms;
  amrex::Vector<int> m_refRatio;
  amrex::Vector<amrex::DistributionMapping> m_dmaps;
  amrex::Vector<amrex::MultiFab> m_data;
};

} // namespace pltfilemanager
} // namespace physics
} // namespace pele
#endif
