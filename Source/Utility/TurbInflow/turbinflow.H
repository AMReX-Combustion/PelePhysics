#ifndef TURBINFLOW_H
#define TURBINFLOW_H

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

namespace pele::physics::turbinflow {

struct TurbParm
{
  // Turbulence data file
  std::string m_turb_file;

  // Verbosity
  int verbose = 1;

  // Geometry of the turbulent data
  amrex::GpuArray<int, 3> npboxcells = {{0}};
  amrex::GpuArray<amrex::Real, 3> pboxlo = {{0.0}};
  amrex::GpuArray<amrex::Real, 3> dx = {{0.0}};
  amrex::GpuArray<amrex::Real, 3> dxinv = {{0.0}};
  amrex::GpuArray<amrex::Real, 3> pboxsize = {{0.0}};
  int dir;                       // Direction where to use this TurbParm
  amrex::Orientation::Side side; // Side on which to use this TurbParm
  amrex::Real time_shift =
    0.0; // Enable to offset the location in the turb file

  int nplane = 32; // Number of turb planes stored at once in memory
  amrex::FArrayBox* sdata = nullptr; // Currently loaded chunk of turb data
  amrex::Real szlo = -1.0e12;        // Position of the current chunk low plane
  amrex::Real szhi = -1.0e11;        // Position of the current chunk high plane

  bool isswirltype = false; // Unused: enable rotating velocity
  amrex::Real turb_scale_loc =
    1.; // Spatial scale relating turbfile geom -> actual case scale
  amrex::Real turb_scale_vel = 1.; // Velocity scaling
  amrex::Real turb_conv_vel =
    1.; // Mean inlet velocity, used to move through data planes

  amrex::Vector<long> m_offset;
  int kmax;     // Number of plane in Turbfile
  long* offset; // Binary offset of each data plane in Turbfile
  long offset_size;
};

struct TurbInflow
{
public:
  TurbInflow() = default;

  ~TurbInflow() = default;

  void init(amrex::Geometry const& geom);

  void add_turb(
    amrex::Box const& bx,
    amrex::FArrayBox& data,
    const int dcomp,
    amrex::Geometry const& geom,
    const amrex::Real time,
    const int dir,
    const amrex::Orientation::Side& side);

  bool is_initialized() const { return turbinflow_initialized; }

  static void read_turb_planes(TurbParm& a_tp, amrex::Real z);

  static void read_one_turb_plane(TurbParm& a_tp, int iplane, int k);

  static void fill_turb_plane(
    TurbParm& a_tp,
    const amrex::Vector<amrex::Real>& x,
    const amrex::Vector<amrex::Real>& y,
    amrex::Real z,
    amrex::FArrayBox& v);

  static void set_turb(
    int normDir,
    int transDir1,
    int transDir2,
    amrex::FArrayBox& v,
    amrex::FArrayBox& data,
    const int dcomp);

private:
  amrex::Vector<TurbParm> tp;
  bool turbinflow_initialized = false;
};
} // namespace pele::physics::turbinflow

#endif
