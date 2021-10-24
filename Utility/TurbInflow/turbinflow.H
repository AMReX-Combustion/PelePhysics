#ifndef _turbinflow_H_
#define _turbinflow_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Geometry.H>

struct TurbParmHost
{
  amrex::Gpu::DeviceVector<long> offset_dv;
  std::string turb_file = "";
};

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
  bool turbinflow_initialized = false;
  bool turbinflow_planes_initialized = false;
  int kmax;
  long* offset;
  long offset_size;
  TurbParmHost* tph = nullptr;
};

void init_turbinflow(
  const std::string& turb_file,
  amrex::Real turb_scale_loc,
  amrex::Real turb_scale_vel,
  const amrex::Vector<amrex::Real>& turb_center,
  amrex::Real turb_conv_vel,
  int nplane,
  TurbParm& tp);

void add_turb(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const int dir,
  const amrex::Orientation::Side& side,
  TurbParm& tp);

void fill_with_turb(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  amrex::Geometry const& geom,
  TurbParm& tp);

void set_turb(amrex::FArrayBox& v,
              amrex::FArrayBox& data);
#endif
