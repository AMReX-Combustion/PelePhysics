#include "Transport.H"

void
transport_init()
{
  transport_params::init();
}

void
transport_close()
{
  transport_params::finalize();
}

AMREX_GPU_DEVICE
void
get_transport_coeffs(
  amrex::Box const& bx,
  amrex::Array4<const amrex::Real> const& Y_in,
  amrex::Array4<const amrex::Real> const& T_in,
  amrex::Array4<const amrex::Real> const& Rho_in,
  amrex::Array4<amrex::Real> const& D_out,
  amrex::Array4<amrex::Real> const& mu_out,
  amrex::Array4<amrex::Real> const& xi_out,
  amrex::Array4<amrex::Real> const& lam_out)
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
          Ddiag, muloc, xiloc, lamloc);

        //   mu, xi and lambda are stored after D in the diffusion multifab
        for (int n = 0; n < NUM_SPECIES; ++n) {
          D_out(i, j, k, n) = Ddiag[n];
        }

        mu_out(i, j, k)  = muloc;
        xi_out(i, j, k)  = xiloc;
        lam_out(i, j, k) = lamloc;
      }
    }
  }
}
