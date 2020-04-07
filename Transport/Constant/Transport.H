#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX.H>
#include <AMReX_Gpu.H>

#include "mechanism.h"
#include "chemistry_file.H"
#include "IndexDefines.H"
#include "TransportParams.H"

void transport_init();
void transport_close();

AMREX_GPU_DEVICE
void get_transport_coeffs(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<amrex::Real> const& D);

using namespace transport_params;

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
transport(
  bool wtr_get_xi,
  bool wtr_get_mu,
  bool wtr_get_lam,
  bool wtr_get_Ddiag,
  amrex::Real& Tloc,
  amrex::Real& rholoc,
  amrex::Real* Yloc,
  amrex::Real* Ddiag,
  amrex::Real& mu,
  amrex::Real& xi,
  amrex::Real& lam)

{

  if (wtr_get_lam) {
    lam = const_conductivity;
  }

  if (wtr_get_mu) {
    mu = const_viscosity;
  }

  if (wtr_get_xi) {
    xi = const_bulk_viscosity;
  }

  if (wtr_get_Ddiag) {
    for (int i = 0; i < NUM_SPECIES; ++i) {
      Ddiag[i] = const_diffusivity;
    }
  }
}

#endif
