#ifndef _MAIN_H
#define _MAIN_H
#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_FArrayBox.H>

#include "PltFileManager.H"

AMREX_ENUM(InputType, turb_box, diag_frame_planes, periodic_plt);

void
turbinflow_from_diag_frame_planes(std::ofstream& ifsd, std::ofstream& ifsh);

void turbinflow_from_periodic_plt(std::ofstream& ifsd, std::ofstream& ifsh);

void turbinflow_from_turb_box(std::ofstream& ifsd, std::ofstream& ifsh);

amrex::IntVect index_mapper(
  const amrex::IntVect& idx_in,
  const amrex::IntVect& dim_map,
  const amrex::Box& domain,
  const amrex::IntVect& periodic);

void read_binary(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<amrex::Real>& data);

void read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data);

struct HITData
{
  int input_ncell = 0;
  amrex::Real urms0 = 1.0;
  amrex::Real uin_norm = 1.0;
  amrex::Real Linput = 0.0;
  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;
};

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fillVelFab(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& vfab,
  amrex::GeometryData const& /*geomdata*/,
  HITData const& a_data)
{
  int inSize = a_data.input_ncell;
  //
  // Fill Velocity
  //
  vfab(i, j, k, 0) = a_data.d_uinput[i + inSize * (j + inSize * k)];
  vfab(i, j, k, 1) = a_data.d_vinput[i + inSize * (j + inSize * k)];
  vfab(i, j, k, 2) = a_data.d_winput[i + inSize * (j + inSize * k)];
}
#endif
