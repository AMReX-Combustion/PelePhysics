#ifndef HITDATA_H
#define HITDATA_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

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
#endif
