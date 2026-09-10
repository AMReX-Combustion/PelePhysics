#ifndef TURBINFLOW_H
#define TURBINFLOW_H

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

// Feature macro: this PelePhysics provides the TurbInflow::add_turb()
// overload that takes the physical positions of the inflow-face cell
// centres explicitly, rather than deriving them from an amrex::Geometry.
//
// A solver whose AMReX grid is a uniform computational (Xi) grid carrying a
// coordinate mapping -- e.g. PeleLMeX's geometry.mesh_mapping -- cannot use
// the Geometry-based overload: the turbulence file is indexed by PHYSICAL
// position, while ProbLo + (i+0.5)*CellSize is the Xi position.  Consumers
// #ifdef on this macro so that the PelePhysics and solver-side changes can
// land independently.
#define PELEPHYSICS_TURBINFLOW_HAS_COORD_ADDTURB 1

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
  bool istimeplanes = false; // false indicates turb data is periodic normal to
                             // flow dir, true requires time array in turb data

  bool extrap_nonperiodic =
    false; // if true, allow interpolation where stencil touches ghost cells for
           // nonperiodic directions
  bool tile_periodic = false; // if true, turbulence is patched over entire face
                              // with a periodic tiled pattern
  bool time_periodic = false; // if true, turbulence data is looped in time
  TurbInterpType interp_type =
    TurbInterpType::quadratic; // linear required for non-periodic directions

  amrex::Real turb_scale_loc =
    1.; // Spatial scale relating turbfile geom -> actual case scale
  amrex::Real turb_scale_vel = 1.; // Velocity scaling
  amrex::Real turb_conv_vel =
    1.; // Mean inlet velocity, used to move through data planes

  amrex::Real convected_distance =
    -1.; // If >= 0, used directly as the through-plane sampling position
         // (physical convected distance, e.g. the time integral of the inlet
         // velocity), bypassing turb_conv_vel * time. Negative => disabled.

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

  //! Uniform-grid entry point: the physical position of a cell centre is
  //! ProbLo + (i+0.5)*CellSize.  Builds those coordinates and forwards to
  //! the overload below, which holds the single implementation.
  void add_turb(
    amrex::Box const& bx,
    amrex::FArrayBox& data,
    const int dcomp,
    amrex::Geometry const& geom,
    const amrex::Real time,
    const int dir,
    const amrex::Orientation::Side& side);

  /**
   * \brief Add turbulence to \p bx using caller-supplied physical
   *        coordinates for the two transverse directions.
   *
   * For a solver on a non-uniform or mapped grid this is the correct entry
   * point: the turbulence file is sampled by physical position, which the
   * caller (who owns the mapping) is the only one able to compute.
   *
   * \p x_phys and \p y_phys are the physical cell-centre coordinates along
   * the transverse directions \c tdir1 and \c tdir2 given by
   * transverseDirs(dir, ...), indexed from \c bx.smallEnd(tdir1) and
   * \c bx.smallEnd(tdir2) respectively, and must span the corresponding
   * extents of \p bx.  They are in case units; the per-TurbParm
   * \c turb_scale_loc is applied internally, exactly as for the
   * Geometry-based overload.
   *
   * \p domain is the (Xi-space) problem domain Box, used only to locate the
   * ghost plane being filled.
   */
  void add_turb(
    amrex::Box const& bx,
    amrex::FArrayBox& data,
    const int dcomp,
    amrex::Box const& domain,
    const amrex::Vector<amrex::Real>& x_phys,
    const amrex::Vector<amrex::Real>& y_phys,
    const amrex::Real time,
    const int dir,
    const amrex::Orientation::Side& side);

  /**
   * \brief The two transverse directions associated with injection
   *        direction \p dir, in the order the turbulence data uses them.
   *
   * Exposed so that callers building coordinate vectors for the overload
   * above cannot disagree with the internal convention.
   */
  static AMREX_FORCE_INLINE void
  transverseDirs(const int dir, int& tdir1, int& tdir2) noexcept
  {
    tdir1 = (dir != 0) ? 0 : 1;
    tdir2 = (dir != 0) ? ((dir == 2) ? 1 : 2) : 2;
  }

  /**
   * \brief Transverse cell spacings of the turbulence file feeding the
   *        (\p dir, \p side) face, expressed in case units.
   *
   * Returns false when no turb file acts on that face.  Intended for
   * resolution-consistency diagnostics: the ratio of the target grid's
   * local spacing to this value is what decides whether the injected
   * field can fill the grid's resolved scales.
   */
  bool file_transverse_dx(
    const int dir,
    const amrex::Orientation::Side& side,
    amrex::Real& dx_tdir1,
    amrex::Real& dx_tdir2) const;

  bool is_initialized() const { return turbinflow_initialized; }

  // Override the through-plane sampling position with a supplied physical
  // convected distance (e.g. the time-integrated inlet velocity). Applies to
  // all loaded turb planes. Pass a negative value to revert to the default
  // turb_conv_vel * time behavior.
  void set_convected_distance(amrex::Real dist)
  {
    for (auto& t : tp) {
      t.convected_distance = dist;
    }
  }

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
