#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include "chemistry_file.H"
#include "TransportParams.H"

namespace pele {
namespace physics {
namespace transport {

template <typename EOSType>
struct InitTransport
{
  void operator()()
  {
    TransParm trans_parm;

    // Default
    trans_parm.Prandtl_number = 0.7;
    trans_parm.viscosity_mu_ref = 17.16;
    trans_parm.viscosity_T_ref = 273.15;
    trans_parm.viscosity_S = 110.4;
    trans_parm.const_bulk_viscosity = 0.0;
    trans_parm.const_diffusivity = 1.0;

    // User-specified
    amrex::ParmParse pp("transport");
    pp.query("Prandtl_number", trans_parm.Prandtl_number);
    pp.query("viscosity_mu_ref", trans_parm.viscosity_mu_ref);
    pp.query("viscosity_T_ref", trans_parm.viscosity_T_ref);
    pp.query("viscosity_S", trans_parm.viscosity_S);
    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
    pp.query("const_diffusivity", trans_parm.const_diffusivity);

    /* GPU */
    trans_parm_g =
      (TransParm*)amrex::The_Device_Arena()->alloc(sizeof(trans_parm));
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
    std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
  }
};

template <typename EOSType>
struct CloseTransport
{
  void operator()() { amrex::The_Device_Arena()->free(trans_parm_g); }
};

struct Sutherland
{
  using transport_type = Sutherland;

  static std::string identifier() { return "Sutherland"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void transport(
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
    amrex::Real& lam,
    TransParm const* trans_parm)
  {

    if (wtr_get_mu) {
      mu = trans_parm->viscosity_mu_ref *
           std::pow((Tloc / trans_parm->viscosity_T_ref), 1.5) *
           (trans_parm->viscosity_T_ref + trans_parm->viscosity_S) /
           (Tloc + trans_parm->viscosity_S);
    }

    if (wtr_get_lam) {
      amrex::Real muloc =
        trans_parm->viscosity_mu_ref *
        std::pow((Tloc / trans_parm->viscosity_T_ref), 1.5) *
        (trans_parm->viscosity_T_ref + trans_parm->viscosity_S) /
        (Tloc + trans_parm->viscosity_S);

      amrex::Real Cpmix;
      CKCPBS(&Tloc, Yloc, &Cpmix);
      lam = muloc * Cpmix / trans_parm->Prandtl_number;
    }

    if (wtr_get_xi) {
      xi = trans_parm->const_bulk_viscosity;
    }

    if (wtr_get_Ddiag) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        Ddiag[i] = trans_parm->const_diffusivity;
      }
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_transport_coeffs(
    amrex::Box const& bx,
    amrex::Array4<const amrex::Real> const& Y_in,
    amrex::Array4<const amrex::Real> const& T_in,
    amrex::Array4<const amrex::Real> const& Rho_in,
    amrex::Array4<amrex::Real> const& D_out,
    amrex::Array4<amrex::Real> const& mu_out,
    amrex::Array4<amrex::Real> const& xi_out,
    amrex::Array4<amrex::Real> const& lam_out,
    TransParm const* trans_parm)
  {
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    bool wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag;

    wtr_get_xi = true;
    wtr_get_mu = true;
    wtr_get_lam = true;
    wtr_get_Ddiag = true;

    amrex::Real T;
    amrex::Real rho;
    amrex::Real massloc[NUM_SPECIES];

    amrex::Real muloc, xiloc, lamloc;
    amrex::Real Ddiag[NUM_SPECIES];

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {

          T = T_in(i, j, k);
          rho = Rho_in(i, j, k);
          for (int n = 0; n < NUM_SPECIES; ++n) {
            massloc[n] = Y_in(i, j, k, n);
          }

          transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, T, rho, massloc,
            Ddiag, muloc, xiloc, lamloc, trans_parm);

          // mu, xi and lambda are stored after D in the diffusion multifab
          for (int n = 0; n < NUM_SPECIES; ++n) {
            D_out(i, j, k, n) = Ddiag[n];
          }

          mu_out(i, j, k) = muloc;
          xi_out(i, j, k) = xiloc;
          lam_out(i, j, k) = lamloc;
        }
      }
    }
  }
  template <class... Args>
  AMREX_GPU_HOST_DEVICE Sutherland(Args...)
  {
  }
};
} // namespace transport
using TransportType = transport::Sutherland;
} // namespace physics
} // namespace pele
#endif
