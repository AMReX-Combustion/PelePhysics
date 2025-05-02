#ifndef MAIN_K_H
#define MAIN_K_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <HITData.H>
#include <Utilities.H>

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
