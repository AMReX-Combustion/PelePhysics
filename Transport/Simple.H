#ifndef _SIMPLETRANSPORT_H_
#define _SIMPLETRANSPORT_H_

#include "TransportParams.H"
#include "PhysicsConstants.H"

namespace pele {
namespace physics {
namespace transport {

template <typename EOSType>
struct NonIdealChungCorrections
{
  template <class... Args>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void operator()(Args... /*unused*/)
  {
  }
};

template <>
struct NonIdealChungCorrections<eos::SRK>
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const bool wtr_get_mu,
    const bool wtr_get_lam,
    const amrex::Real Temp,
    const amrex::Real* Xloc,
    const amrex::Real rholoc,
    amrex::Real& wbar,
    amrex::Real& mu,
    amrex::Real& lam,
    TransParm<eos::SRK, SimpleTransport> const* trans_parm)
  {
    amrex::Real sigma_M_3 = 0.0;
    amrex::Real Epsilon_M = 0.0;
    amrex::Real Omega_M = 0.0;
    amrex::Real MW_m = 0.0;
    amrex::Real DP_m_4 = 0.0;
    amrex::Real KappaM = 0.0;

    // Note: all the square roots are precalculated for efficiency

    for (int i = 0; i < NUM_SPECIES; ++i) {
      for (int j = 0; j < NUM_SPECIES; ++j) {
        const amrex::Real Xij = Xloc[i] * Xloc[j];
        const int idx = i * NUM_SPECIES + j;

        const amrex::Real sqrtT2 = trans_parm->sqrtT2ij[idx];
        const amrex::Real T2 = sqrtT2 * sqrtT2;
        const amrex::Real T3 = T2 * sqrtT2;
        sigma_M_3 += Xij * T3;

        const amrex::Real Epsilon_ij = trans_parm->sqrtEpsilonij[idx];
        Epsilon_M += Xij * Epsilon_ij * T3;

        const amrex::Real Omega_ij =
          0.5 * (trans_parm->omega[i] + trans_parm->omega[j]);
        Omega_M += Xij * Omega_ij * T3;

        MW_m += Xij * Epsilon_ij * T2 * trans_parm->sqrtMWij[idx];

        DP_m_4 += Xij * trans_parm->dip[i] * trans_parm->dip[i] *
                  trans_parm->dip[j] * trans_parm->dip[j] / (T3 * Epsilon_ij);

        KappaM += Xij * trans_parm->sqrtKappaij[idx];
      }
    }

    MW_m *= MW_m;
    const amrex::Real sigma_M = std::cbrt(sigma_M_3);
    const amrex::Real InvSigma1 = 1.0 / sigma_M;
    const amrex::Real InvSigma3 = 1.0 / sigma_M_3;
    Epsilon_M *= InvSigma3;
    const amrex::Real Tstar = Temp / Epsilon_M;
    const amrex::Real Tcm = 1.2593 * Epsilon_M;
    const amrex::Real Vcm = 1.8887 * sigma_M_3;
    Omega_M *= InvSigma3;
    MW_m = MW_m * InvSigma3 * InvSigma1 / (Epsilon_M * Epsilon_M);
    DP_m_4 = DP_m_4 * sigma_M_3 * Epsilon_M;
    const amrex::Real DP_red_4 =
      297.2069113e6 * DP_m_4 / (Vcm * Vcm * Tcm * Tcm);
    const amrex::Real y = Vcm * rholoc / (6.0 * wbar);
    amrex::Real A[10];
    for (int i = 0; i < 10; ++i) {
      A[i] = trans_parm->Afac[4 * i] + trans_parm->Afac[4 * i + 1] * Omega_M +
             trans_parm->Afac[4 * i + 2] * DP_red_4 +
             trans_parm->Afac[4 * i + 3] * KappaM;
    }
    const amrex::Real G1 =
      (1.0 - 0.5 * y) / ((1.0 - y) * (1.0 - y) * (1.0 - y));
    const amrex::Real G2 = (A[0] * (1.0 - std::exp(-A[3] * y)) / y +
                            A[1] * G1 * std::exp(A[4] * y) + A[2] * G1) /
                           (A[0] * A[3] + A[1] + A[2]);
    const amrex::Real eta_P =
      (36.344e-6 * std::sqrt(MW_m * Tcm) / std::pow(Vcm, 2.0 / 3.0)) * A[6] *
      y * y * G2 * exp(A[7] + A[8] / Tstar + A[9] / (Tstar * Tstar));

    // Set nonideal viscosity
    if (wtr_get_mu) {
      mu = mu * (1.0 / G2 + A[5] * y) + eta_P;
    }

    amrex::Real B[7];
    for (int i = 0; i < 7; ++i) {
      B[i] = trans_parm->Bfac[i * 4] + trans_parm->Bfac[i * 4 + 1] * Omega_M +
             trans_parm->Bfac[i * 4 + 2] * DP_red_4 +
             trans_parm->Bfac[i * 4 + 3] * KappaM;
    }
    amrex::Real H2 = (B[0] * (1.0 - std::exp(-B[3] * y)) / y +
                      B[1] * G1 * std::exp(B[4] * y) + B[2] * G1) /
                     (B[0] * B[3] + B[1] + B[2]);

    amrex::Real lambda_p = 3.039e-4 * std::sqrt(Tcm / MW_m) /
                           std::pow(Vcm, 2.0 / 3.0) * B[6] * y * y * H2 *
                           std::sqrt(Tstar);

    lambda_p *= 4.184e+7; // erg/(cm s K)
    const amrex::Real beta = 1.0 / H2 + B[5] * y;

    // Set nonideal conductivity
    if (wtr_get_lam) {
      lam = lam * beta + lambda_p;
    }
  }
};

template <typename EOSType>
struct BinaryDiff
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const amrex::Real* Xloc,
    const amrex::Real* Yloc,
    const amrex::Real* logT,
    const amrex::Real /*rholoc*/,
    const amrex::Real Tloc,
    amrex::Real* Ddiag,
    TransParm<EOSType, SimpleTransport> const* tparm)
  {
    const amrex::Real scale = Constants::PATM / (Constants::RU * Tloc);
    for (int i = 0; i < NUM_SPECIES; ++i) {
      amrex::Real term1 = 0.0;
      amrex::Real term2 = 0.0;
      for (int j = 0; j < NUM_SPECIES; ++j) {
        if (i != j) {
          const int four_idx_ij = 4 * (i + NUM_SPECIES * j);
          const amrex::Real dbintemp =
            tparm->fitdbin[four_idx_ij] +
            tparm->fitdbin[1 + four_idx_ij] * logT[0] +
            tparm->fitdbin[2 + four_idx_ij] * logT[1] +
            tparm->fitdbin[3 + four_idx_ij] * logT[2];
          term1 += Yloc[j];
          term2 += Xloc[j] * std::exp(-dbintemp);
        }
      }
      Ddiag[i] = tparm->wt[i] * term1 / term2 * scale;
    }
  }
};

template <>
struct BinaryDiff<eos::SRK>
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const amrex::Real* Xloc,
    const amrex::Real* Yloc,
    const amrex::Real* logT,
    const amrex::Real rholoc,
    const amrex::Real Tloc,
    amrex::Real* Ddiag,
    TransParm<eos::SRK, SimpleTransport> const* tparm)
  {
    for (int i = 0; i < NUM_SPECIES; ++i) {
      amrex::Real term1 = 0.0;
      amrex::Real term2 = 0.0;
      for (int j = 0; j < NUM_SPECIES; ++j) {
        if (i != j) {
          const int idx_ij = i + NUM_SPECIES * j;
          amrex::Real dbintemp = tparm->fitdbin[4 * idx_ij] +
                                 tparm->fitdbin[1 + 4 * idx_ij] * logT[0] +
                                 tparm->fitdbin[2 + 4 * idx_ij] * logT[1] +
                                 tparm->fitdbin[3 + 4 * idx_ij] * logT[2];
          dbintemp = std::exp(-dbintemp);

          amrex::Real Upsilonij = 0.0;
          for (int k = 0; k < NUM_SPECIES; ++k) {
            Upsilonij += tparm->Upsilonijk[idx_ij * NUM_SPECIES + k] * Yloc[k];
          }
          Upsilonij = Upsilonij * rholoc * Constants::Avna * M_PI / 12.0 + 1.0;
          dbintemp *= (Constants::RU * Tloc * Upsilonij) / Constants::PATM;
          term1 += Yloc[j];
          term2 += Xloc[j] * dbintemp;
        }
      }
      Ddiag[i] = tparm->wt[i] * term1 / term2;
    }
  }
};

struct SimpleTransport
{
  using transport_type = SimpleTransport;

  static std::string identifier() { return "SimpleTransport"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void comp_pure_bulk(
    amrex::Real Tloc,
    const amrex::Real* muloc,
    amrex::Real* xiloc,
    TransParm<EosType, transport_type> const* tparm)
  {
    amrex::Real cvk[NUM_SPECIES];
    amrex::Real cvkint[NUM_SPECIES];
    amrex::Real cvkrot[NUM_SPECIES];
    amrex::Real FofT[NUM_SPECIES];
    amrex::Real Fnorm[NUM_SPECIES];
    amrex::Real epskoverT, epskoverTstd;
    amrex::Real pi3_2 = std::pow(M_PI, 1.5);

    CKCVMS(&Tloc, cvk);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      if (tparm->nlin[i] == 0) {
        cvkint[i] = 0.0;
        cvkrot[i] = 0.0;
      } else if (tparm->nlin[i] == 1) {
        cvkint[i] = cvk[i] * tparm->wt[i] / Constants::RU - 1.50;
        cvkrot[i] = 1.0;
      } else {
        cvkint[i] = cvk[i] * tparm->wt[i] / Constants::RU - 1.50;
        cvkrot[i] = 1.5;
      }
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      epskoverTstd = tparm->eps[i] / 298.0;
      epskoverT = tparm->eps[i] / Tloc;
      Fnorm[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverTstd) +
                 (2.0 + .50 * M_PI * M_PI) * epskoverTstd +
                 pi3_2 * std::sqrt(epskoverTstd) * epskoverTstd;

      FofT[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverT) +
                (2.0 + .50 * M_PI * M_PI) * epskoverT +
                pi3_2 * std::sqrt(epskoverT) * epskoverT;
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      if (tparm->nlin[i] == 0) {
        xiloc[i] = 0.0;
      } else {
        //   zrot/crot approximately zint / cint by assuming vibrational
        //   internal energy is small cvkrot is scaled by wk / Ru = mk / kb
        //   relative to standard specific cv
        xiloc[i] =
          0.250 * M_PI *
          (cvkint[i] / (cvkint[i] + 1.50) * cvkint[i] / (cvkint[i] + 1.50)) *
          tparm->zrot[i] / cvkrot[i] * Fnorm[i] / FofT[i] * muloc[i];
      }
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void transport(
    const bool wtr_get_xi,
    const bool wtr_get_mu,
    const bool wtr_get_lam,
    const bool wtr_get_Ddiag,
    const amrex::Real Tloc,
    const amrex::Real rholoc,
    amrex::Real* Yloc,
    amrex::Real* Ddiag,
    amrex::Real& mu,
    amrex::Real& xi,
    amrex::Real& lam,
    TransParm<EosType, transport_type> const* tparm)
  {
    amrex::Real trace = 1.e-15;
    amrex::Real Xloc[NUM_SPECIES];
    amrex::Real muloc[NUM_SPECIES];
    amrex::Real xiloc[NUM_SPECIES];
    amrex::Real lamloc[NUM_SPECIES];
    amrex::Real logT[NUM_FIT - 1];

    logT[0] = std::log(Tloc);
    logT[1] = logT[0] * logT[0];
    logT[2] = logT[0] * logT[1];

    amrex::Real sum = 0.;

    for (int i = 0; i < NUM_SPECIES; ++i) {
      sum += Yloc[i];
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      Yloc[i] += trace * (sum / NUM_SPECIES - Yloc[i]);
    }

    amrex::Real wbar = 0.;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      wbar += Yloc[i] * tparm->iwt[i];
    }
    wbar = 1.0 / wbar;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      Xloc[i] = Yloc[i] * wbar * tparm->iwt[i];
    }
    if (wtr_get_mu) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        muloc[i] = tparm->fitmu[4 * i] + tparm->fitmu[1 + 4 * i] * logT[0] +
                   tparm->fitmu[2 + 4 * i] * logT[1] +
                   tparm->fitmu[3 + 4 * i] * logT[2];
        muloc[i] = std::exp(muloc[i]);
      }

      mu = 0.0;

      for (int i = 0; i < NUM_SPECIES; ++i) {
        mu += Xloc[i] *
              (muloc[i] * muloc[i] * muloc[i] * muloc[i] * muloc[i] * muloc[i]);
      }
      mu = std::cbrt(std::sqrt(mu)); // mu = (sum[ Xloc_i * muloc_i^6])^(1/6)

      //  assumption that we only get bulk viscosity if we are already getting
      //  shear viscosity

      if (wtr_get_xi) {
        comp_pure_bulk(Tloc, muloc, xiloc, tparm);
        xi = 0.0;
        for (int i = 0; i < NUM_SPECIES; ++i) {
          xi += Xloc[i] * (std::sqrt(std::sqrt(
                            xiloc[i] * xiloc[i] *
                            xiloc[i]))); // xi = sum[Xloc_i*xiloc_i^(3/4)];
        }
        xi = std::cbrt(
          xi * xi * xi * xi); // xi = (sum[Xloc_i*xiloc_i^(3/4)])^(4/3)
      }
    }

    if (wtr_get_lam) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        lamloc[i] = tparm->fitlam[4 * i] + tparm->fitlam[1 + 4 * i] * logT[0] +
                    tparm->fitlam[2 + 4 * i] * logT[1] +
                    tparm->fitlam[3 + 4 * i] * logT[2];
        lamloc[i] = std::exp(lamloc[i]);
      }

      lam = 0.;
      for (int i = 0; i < NUM_SPECIES; ++i) {
        lam += Xloc[i] * std::sqrt(std::sqrt(
                           lamloc[i])); // lam = sum[Xloc_i * (lamloc_i)^(1/4)]
      }
      lam = lam * lam * lam * lam; // lam = lam^4
    }
    // Add nonideal corrections if necessary
    NonIdealChungCorrections<EosType>()(
      wtr_get_mu, wtr_get_lam, Tloc, Xloc, rholoc, wbar, mu, lam, tparm);

    if (wtr_get_Ddiag) {
      BinaryDiff<EosType>()(Xloc, Yloc, logT, rholoc, Tloc, Ddiag, tparm);
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
    TransParm<EosType, transport_type> const* tparm)
  {
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const bool wtr_get_xi = true;
    const bool wtr_get_mu = true;
    const bool wtr_get_lam = true;
    const bool wtr_get_Ddiag = true;

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {

          const amrex::Real T = T_in(i, j, k);
          const amrex::Real rho = Rho_in(i, j, k);
          amrex::Real massloc[NUM_SPECIES] = {0.0};
          for (int n = 0; n < NUM_SPECIES; ++n) {
            massloc[n] = Y_in(i, j, k, n);
          }

          amrex::Real muloc, xiloc, lamloc;
          amrex::Real Ddiag[NUM_SPECIES] = {0.0};
          transport(
            wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, T, rho, massloc,
            Ddiag, muloc, xiloc, lamloc, tparm);

          //   mu, xi and lambda are stored after D in the diffusion multifab
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
  AMREX_GPU_HOST_DEVICE SimpleTransport(Args... /*unused*/)
  {
  }
};
} // namespace transport
} // namespace physics
} // namespace pele
#endif
