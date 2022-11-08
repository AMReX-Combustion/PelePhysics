#ifndef CONSTTRANSPORT_H
#define CONSTTRANSPORT_H

#include "TransportParams.H"

namespace pele {
namespace physics {
namespace transport {
struct ConstTransport
{
  using transport_type = ConstTransport;

  static std::string identifier() { return "ConstTransport"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void transport(
    const bool wtr_get_xi,
    const bool wtr_get_mu,
    const bool wtr_get_lam,
    const bool wtr_get_Ddiag,
    const bool wtr_get_chi,
    const amrex::Real /*Tloc*/,
    const amrex::Real /*rholoc*/,
    amrex::Real* /*Yloc*/,
    amrex::Real* Ddiag,
    amrex::Real* chi_mix,
    amrex::Real& mu,
    amrex::Real& xi,
    amrex::Real& lam,
    TransParm<EosType, transport_type> const* tparm)
  {

    if (wtr_get_lam) {
      lam = tparm->const_conductivity;
    }

    if (wtr_get_mu) {
      mu = tparm->const_viscosity;
    }

    if (wtr_get_xi) {
      xi = tparm->const_bulk_viscosity;
    }

    if (wtr_get_Ddiag) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        Ddiag[i] = tparm->const_diffusivity;
      }
    }
    if (wtr_get_chi) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        chi_mix[i] = tparm->const_thermal_diffusion_ratio;
      }
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE static void get_transport_coeffs(
    amrex::Box const& bx,
    amrex::Array4<const amrex::Real> const& /*Y_in*/,
    amrex::Array4<const amrex::Real> const& /*T_in*/,
    amrex::Array4<const amrex::Real> const& /*Rho_in*/,
    amrex::Array4<amrex::Real> const& D_out,
    amrex::Array4<amrex::Real> const& chi_out,
    amrex::Array4<amrex::Real> const& mu_out,
    amrex::Array4<amrex::Real> const& xi_out,
    amrex::Array4<amrex::Real> const& lam_out,
    TransParm<EosType, transport_type> const* tparm)
  {
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const bool wtr_get_xi = true;
    const bool wtr_get_mu = true;
    const bool wtr_get_lam = true;
    const bool wtr_get_Ddiag = true;
    const bool wtr_get_chi = tparm->use_soret;

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {

          amrex::Real muloc, xiloc, lamloc;
          amrex::Real Ddiag[NUM_SPECIES] = {0.0};
          amrex::Real chi_mix[NUM_SPECIES] = {0.0};
          transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi,
            0.0, 0.0, nullptr, Ddiag, chi_mix, muloc, xiloc, lamloc, tparm);

          // mu, xi and lambda are stored after D in the diffusion multifab
          for (int n = 0; n < NUM_SPECIES; ++n) {
            D_out(i, j, k, n) = Ddiag[n];
            chi_out(i, j, k, n) = chi_mix[n];
          }

          mu_out(i, j, k) = muloc;
          xi_out(i, j, k) = xiloc;
          lam_out(i, j, k) = lamloc;
        }
      }
    }
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE explicit ConstTransport(Args... /*unused*/)
  {
  }
};
} // namespace transport
} // namespace physics
} // namespace pele
#endif
