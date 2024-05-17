#ifndef _MANIFOLDTRANSPORT_H_
#define _MANIFOLDTRANSPORT_H_

#include "TransportParams.H"

namespace pele::physics::transport {

struct ManifoldTransport
{
  using transport_type = ManifoldTransport;

  static std::string identifier() { return "Manifold"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void transport(
    const bool wtr_get_xi,
    const bool wtr_get_mu,
    const bool wtr_get_lam,
    const bool wtr_get_Ddiag,
    const bool wtr_get_chi,
    const amrex::Real /*Tloc*/,
    const amrex::Real /*rholoc*/,
    amrex::Real* Yloc,
    amrex::Real* Ddiag,
    amrex::Real* /*chi_mix*/,
    amrex::Real& mu,
    amrex::Real& /*xi*/,
    amrex::Real& /*lam*/,
    TransParm<EosType, transport_type> const* tparm)
  {
    pele::physics::ManifoldFunc* manfunc;
    if (tparm->manf_data->manmodel == ManifoldModel::TABLE) {
      const TabFuncData* tf_data =
        static_cast<const TabFuncData*>(tparm->manf_data);
      manfunc = new pele::physics::TabFunc(tf_data);
    } else {
      const NNFuncData* nnf_data =
        static_cast<const NNFuncData*>(tparm->manf_data);
      manfunc = new pele::physics::NNFunc(nnf_data);
    }

    if (wtr_get_mu) {
      manfunc->get_value(tparm->idx_mu, Yloc, mu);
    }

    if (wtr_get_lam) {
      // could be added to table if desired
      amrex::Abort(
        "Use of thermal conductivity is not supported for ManifoldTransport");
    }

    if (wtr_get_xi) {
      // could be added to table if desired
      amrex::Abort(
        "Use of bulk viscosity is not supported for ManifoldTransport");
    }

    if (wtr_get_chi) {
      // could be added to table if desired
      amrex::Abort(
        "Use of soret diffusion is not supported for ManifoldTransport");
    }

    if (wtr_get_Ddiag) {
      amrex::Real Dtemp{0.0};
      manfunc->get_value(tparm->idx_rhoD, Yloc, Dtemp);
      for (int i = 0; i < NUM_SPECIES - 1; ++i) {
        Ddiag[i] = Dtemp;
      }
      Ddiag[NUM_SPECIES - 1] = 0.0; // Placeholder for density
    }

    delete manfunc;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void get_transport_coeffs(
    amrex::Box const& bx,
    amrex::Array4<const amrex::Real> const& Y_in,
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

    const bool wtr_get_xi = false;
    const bool wtr_get_mu = true;
    const bool wtr_get_lam = false;
    const bool wtr_get_Ddiag = true;
    const bool wtr_get_chi = false;

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {

          amrex::Real massloc[NUM_SPECIES] = {0.0};
          for (int n = 0; n < NUM_SPECIES; ++n) {
            massloc[n] = Y_in(i, j, k, n);
          }

          amrex::Real muloc = 0.0, xiloc = 0.0, lamloc = 0.0;
          amrex::Real Ddiag[NUM_SPECIES] = {0.0};
          amrex::Real chi_loc[NUM_SPECIES] = {0.0};
          transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi,
            0.0, 0.0, massloc, Ddiag, chi_loc, muloc, xiloc, lamloc, tparm);

          // mu, xi and lambda are stored after D in the diffusion multifab
          for (int n = 0; n < NUM_SPECIES; ++n) {
            D_out(i, j, k, n) = Ddiag[n];
            chi_out(i, j, k, n) = chi_loc[n];
          }

          mu_out(i, j, k) = muloc;
          xi_out(i, j, k) = xiloc;
          lam_out(i, j, k) = lamloc;
        }
      }
    }
  }
  template <class... Args>
  AMREX_GPU_HOST_DEVICE ManifoldTransport(Args...)
  {
  }
};

} // namespace pele::physics::transport
#endif
