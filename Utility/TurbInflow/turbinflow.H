#ifndef _turbinflow_H_
#define _turbinflow_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

namespace pele {
namespace physics {
namespace turbinflow {

struct TurbParm
{
  amrex::GpuArray<int, AMREX_SPACEDIM> npboxcells = {{0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pboxlo = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv = {{0.0}};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pboxsize = {{0.0}};
  int nplane = 32;
  amrex::FArrayBox* sdata = nullptr;
  amrex::Real szlo = 0.0;
  amrex::Real szhi = 0.0;
  bool isswirltype = false;
  amrex::Real turb_scale_loc = 1.;
  amrex::Real turb_scale_vel = 1.;
  amrex::Real turb_conv_vel = 1.;
  int kmax;
  long* offset;
  long offset_size;
};

struct TurbInflow
{
public:
  TurbInflow() = default;
  ;

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

  void read_turb_planes(amrex::Real z);

  void read_one_turb_plane(int iplane, int k);

  void fill_turb_plane(
    const amrex::Vector<amrex::Real>& x,
    const amrex::Vector<amrex::Real>& y,
    amrex::Real z,
    amrex::FArrayBox& v);

  void fill_with_turb(
    amrex::Box const& bx,
    amrex::FArrayBox& data,
    const int dcomp,
    amrex::Geometry const& geom);

  void set_turb(
    int nDir,
    int trDir1,
    int trDir2,
    amrex::FArrayBox& v,
    amrex::FArrayBox& data,
    const int dcomp);

private:
  std::string m_turb_file = "";

  TurbParm tp;

  amrex::Gpu::DeviceVector<long> m_offset_dv;

  bool turbinflow_initialized = false;

  bool turbinflow_planes_initialized = false;
};
} // namespace turbinflow

} // namespace physics
} // namespace pele

#endif
