#ifndef TURBINFLOW_H
#define TURBINFLOW_H

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

namespace pele::physics::turbinflow {

AMREX_ENUM(TurbInterpType, linear, quadratic);

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
  int izlo = -1;                     // Currently loaded low plane
  int izhi = -1;                     // Currently loaded high plane
  amrex::Real szlo = -1.0e12;        // Position of the current chunk low plane
  amrex::Real szhi = -1.0e11;        // Position of the current chunk high plane

  amrex::IntVect periodicity{AMREX_D_DECL(1, 1, 1)};
  bool isswirltype = false; // false indicates turb data is periodic normal to
                            // flow dir, true requires time array in turb data

  bool extrap_nonperiodic =
    false; // if true, allow interpolation where stencil touches ghost cells for
           // nonperiodic directions
  bool tile_periodic = false; // if true, turbulence is patched over entire face
                              // with a periodic tiled pattern
  TurbInterpType interp_type =
    TurbInterpType::quadratic; // linear required for non-periodic directions

  amrex::Real turb_scale_loc =
    1.; // Spatial scale relating turbfile geom -> actual case scale
  amrex::Real turb_scale_vel = 1.; // Velocity scaling
  amrex::Real turb_conv_vel =
    1.; // Mean inlet velocity, used to move through data planes

  amrex::Vector<amrex::Real> planeTimes;

  amrex::Vector<long> offset;
  int kmax; // Number of plane in Turbfile
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
