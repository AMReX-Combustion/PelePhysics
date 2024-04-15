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
  amrex::GeometryData const& geomdata,
  HITData const& a_data)
{

  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  amrex::Real x[3] = {
    prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
    prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
    prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

  // Fill in the velocities
  amrex::Real uinterp[3] = {0.0};

  // Interpolation factors
  amrex::Real mod[3] = {0.0};
  int idx[3] = {0};
  int idxp1[3] = {0};
  amrex::Real slp[3] = {0.0};
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    mod[idim] = std::fmod(x[idim], a_data.Linput);
    locate(a_data.d_xarray, a_data.input_ncell, mod[idim], idx[idim]);
    idxp1[idim] = (idx[idim] + 1) % a_data.input_ncell;
    slp[idim] =
      (mod[idim] - a_data.d_xarray[idx[idim]]) / a_data.d_xdiff[idx[idim]];
  }

  int inSize = a_data.input_ncell;

  const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
  const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
  const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
  const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
  const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
  const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
  const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
  const amrex::Real f7 = slp[0] * slp[1] * slp[2];

  uinterp[0] =
    a_data.d_uinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 +
    a_data.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1 +
    a_data.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2 +
    a_data.d_uinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3 +
    a_data.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4 +
    a_data.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5 +
    a_data.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6 +
    a_data.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

  uinterp[1] =
    a_data.d_vinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 +
    a_data.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1 +
    a_data.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2 +
    a_data.d_vinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3 +
    a_data.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4 +
    a_data.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5 +
    a_data.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6 +
    a_data.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

  uinterp[2] =
    a_data.d_winput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 +
    a_data.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1 +
    a_data.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2 +
    a_data.d_winput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3 +
    a_data.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4 +
    a_data.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5 +
    a_data.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6 +
    a_data.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

  //
  // Fill Velocity
  //
  vfab(i, j, k, 0) = uinterp[0];
  vfab(i, j, k, 1) = uinterp[1];
  vfab(i, j, k, 2) = uinterp[2];
}
#endif
