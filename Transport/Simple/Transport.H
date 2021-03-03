#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <TransportParams.H>

#include "mechanism.h"
#include "chemistry_file.H"
#include "EOS.H"
#include "PhysicsConstants.H"

namespace pele {
namespace physics {
namespace transport {

template <typename EOSType>
struct InitTransport
{
  void operator()()
  {
    trans_parm.trans_wt = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_iwt = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_eps = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_sig = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_dip = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_pol = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_zrot = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);

    trans_parm.trans_fitmu = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
    trans_parm.trans_fitlam = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
    trans_parm.trans_fitdbin = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES * NUM_FIT);
    trans_parm.trans_nlin =
      (int*)amrex::The_Arena()->alloc(sizeof(int) * NUM_SPECIES);

    egtransetWT(trans_parm.trans_wt);
    egtransetEPS(trans_parm.trans_eps);
    egtransetSIG(trans_parm.trans_sig);
    egtransetDIP(trans_parm.trans_dip);
    egtransetPOL(trans_parm.trans_pol);
    egtransetZROT(trans_parm.trans_zrot);
    egtransetNLIN(trans_parm.trans_nlin);
    egtransetCOFETA(trans_parm.trans_fitmu);
    egtransetCOFLAM(trans_parm.trans_fitlam);
    egtransetCOFD(trans_parm.trans_fitdbin);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      trans_parm.trans_iwt[i] = 1. / trans_parm.trans_wt[i];
    }

    trans_parm_g =
      (TransParm*)amrex::The_Device_Arena()->alloc(sizeof(trans_parm));
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
    std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
  }
};

// FIXME remove ifdef when compiling in all the EOS
#ifdef PELEC_USE_SRK
template <>
struct InitTransport<eos::SRK>
{
  void operator()()
  {
    trans_parm.trans_wt = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_iwt = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_eps = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_sig = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_dip = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_pol = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.trans_zrot = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);

    trans_parm.trans_fitmu = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
    trans_parm.trans_fitlam = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_FIT);
    trans_parm.trans_fitdbin = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES * NUM_FIT);
    trans_parm.trans_nlin =
      (int*)amrex::The_Arena()->alloc(sizeof(int) * NUM_SPECIES);

    egtransetWT(trans_parm.trans_wt);
    egtransetEPS(trans_parm.trans_eps);
    egtransetSIG(trans_parm.trans_sig);
    egtransetDIP(trans_parm.trans_dip);
    egtransetPOL(trans_parm.trans_pol);
    egtransetZROT(trans_parm.trans_zrot);
    egtransetNLIN(trans_parm.trans_nlin);
    egtransetCOFETA(trans_parm.trans_fitmu);
    egtransetCOFLAM(trans_parm.trans_fitlam);
    egtransetCOFD(trans_parm.trans_fitdbin);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      trans_parm.trans_iwt[i] = 1. / trans_parm.trans_wt[i];
    }

    // Nonideal transport coefficients computed using Chung's method:
    // Chung, T.H., Ajlan, M., Lee, L.L. and Starling, K.E., 1988.
    // Generalized multiparameter correlation for nonpolar and polar
    // fluid transport properties. Industrial & engineering chemistry
    // research, 27(4), pp.671-679.

    trans_parm.Afac =
      (amrex::Real*)amrex::The_Arena()->alloc(sizeof(amrex::Real) * 10 * 4);
    trans_parm.Bfac =
      (amrex::Real*)amrex::The_Arena()->alloc(sizeof(amrex::Real) * 7 * 4);
    trans_parm.omega = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.Kappai = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES);
    trans_parm.Sigmaij = (amrex::Real*)amrex::The_Arena()->alloc(
      sizeof(amrex::Real) * NUM_SPECIES * NUM_SPECIES);

    // Initialize coefficients of model
    {
      amrex::Real tmp1[NUM_SPECIES], tmp2[NUM_SPECIES], tmp3[NUM_SPECIES];
      GET_CRITPARAMS(tmp1, tmp2, tmp3, trans_parm.omega);
    }

    trans_parm.Afac[0] = 6.32402;
    trans_parm.Afac[1] = 50.41190;
    trans_parm.Afac[2] = -51.68010;
    trans_parm.Afac[3] = 1189.020;
    trans_parm.Afac[4] = 0.12102e-2;
    trans_parm.Afac[5] = -0.11536e-2;
    trans_parm.Afac[6] = -0.62571e-2;
    trans_parm.Afac[7] = 0.37283e-1;
    trans_parm.Afac[8] = 5.28346;
    trans_parm.Afac[9] = 254.209;
    trans_parm.Afac[10] = -168.481;
    trans_parm.Afac[11] = 3898.27;
    trans_parm.Afac[12] = 6.62263;
    trans_parm.Afac[13] = 38.09570;
    trans_parm.Afac[14] = -8.46414;
    trans_parm.Afac[15] = 31.41780;
    trans_parm.Afac[16] = 19.74540;
    trans_parm.Afac[17] = 7.63034;
    trans_parm.Afac[18] = -14.35440;
    trans_parm.Afac[19] = 31.52670;
    trans_parm.Afac[20] = -1.89992;
    trans_parm.Afac[21] = -12.5367;
    trans_parm.Afac[22] = 4.98529;
    trans_parm.Afac[23] = -18.15070;
    trans_parm.Afac[24] = 24.27450;
    trans_parm.Afac[25] = 3.44945;
    trans_parm.Afac[26] = -11.29130;
    trans_parm.Afac[27] = 69.34660;
    trans_parm.Afac[28] = 0.79716;
    trans_parm.Afac[29] = 1.11764;
    trans_parm.Afac[30] = 0.12348e-1;
    trans_parm.Afac[31] = -4.11661;
    trans_parm.Afac[32] = -0.23816;
    trans_parm.Afac[33] = 0.67695e-1;
    trans_parm.Afac[34] = -0.81630;
    trans_parm.Afac[35] = 4.02528;
    trans_parm.Afac[36] = 0.68629e-1;
    trans_parm.Afac[37] = 0.34793;
    trans_parm.Afac[38] = 0.59256;
    trans_parm.Afac[39] = -0.72663;

    trans_parm.Bfac[0] = 2.41657;
    trans_parm.Bfac[1] = 0.74824;
    trans_parm.Bfac[2] = -0.91858;
    trans_parm.Bfac[3] = 121.721;
    trans_parm.Bfac[4] = -0.50924;
    trans_parm.Bfac[5] = -1.50936;
    trans_parm.Bfac[6] = -49.99120;
    trans_parm.Bfac[7] = 69.9834;
    trans_parm.Bfac[8] = 6.61069;
    trans_parm.Bfac[9] = 5.62073;
    trans_parm.Bfac[10] = 64.75990;
    trans_parm.Bfac[11] = 27.0389;
    trans_parm.Bfac[12] = 14.54250;
    trans_parm.Bfac[13] = -8.91387;
    trans_parm.Bfac[14] = -5.63794;
    trans_parm.Bfac[15] = 74.3435;
    trans_parm.Bfac[16] = 0.79274;
    trans_parm.Bfac[17] = 0.82019;
    trans_parm.Bfac[18] = -0.69369;
    trans_parm.Bfac[19] = 6.31734;
    trans_parm.Bfac[20] = -5.86340;
    trans_parm.Bfac[21] = 12.80050;
    trans_parm.Bfac[22] = 9.58926;
    trans_parm.Bfac[23] = -65.5292;
    trans_parm.Bfac[24] = 81.17100;
    trans_parm.Bfac[25] = 114.15800;
    trans_parm.Bfac[26] = -60.84100;
    trans_parm.Bfac[27] = 466.7750;

    for (int i = 0; i < NUM_SPECIES; ++i) {
      for (int j = 0; j < NUM_SPECIES; ++j) {
        trans_parm.Sigmaij[i * NUM_SPECIES + j] =
          0.5 * (trans_parm.trans_sig[i] + trans_parm.trans_sig[j]) *
          1e-8; // converted to cm
      }
    }

    // Initialize Kappa, which has nonzero values only for specific polar
    // species
    for (int i = 0; i < NUM_SPECIES; ++i) {
      trans_parm.Kappai[i] = 0.0;
    }
    {
      amrex::Vector<std::string> spec_names_kappa;
      spec_names_kappa.resize(1);
      spec_names_kappa[0] = "Null";
      CKSYMS_STR(spec_names_kappa);
      for (int i = 0; i < NUM_SPECIES; i++) {
        if (spec_names_kappa[i] == "H2O") {
          trans_parm.Kappai[i] = 0.075908;
          // trans_parm.Kappai[i] = 0.076;
        } else if (spec_names_kappa[i] == "CH3OH") {
          trans_parm.Kappai[i] = 0.215175;
        } else if (spec_names_kappa[i] == "CH3CH2OH") {
          trans_parm.Kappai[i] = 0.174823;
        }
      }
    }

    trans_parm_g =
      (TransParm*)amrex::The_Device_Arena()->alloc(sizeof(trans_parm));
#ifdef AMREX_USE_GPU
    amrex::Gpu::htod_memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#else
    std::memcpy(trans_parm_g, &trans_parm, sizeof(trans_parm));
#endif
  }
};
#endif

template <typename EOSType>
struct CloseTransport
{
  void operator()()
  {
    amrex::The_Device_Arena()->free(trans_parm_g);

    amrex::The_Arena()->free(trans_parm.trans_wt);
    amrex::The_Arena()->free(trans_parm.trans_iwt);
    amrex::The_Arena()->free(trans_parm.trans_eps);
    amrex::The_Arena()->free(trans_parm.trans_sig);
    amrex::The_Arena()->free(trans_parm.trans_dip);
    amrex::The_Arena()->free(trans_parm.trans_pol);
    amrex::The_Arena()->free(trans_parm.trans_zrot);
    amrex::The_Arena()->free(trans_parm.trans_fitmu);
    amrex::The_Arena()->free(trans_parm.trans_fitlam);
    amrex::The_Arena()->free(trans_parm.trans_fitdbin);
    amrex::The_Arena()->free(trans_parm.trans_nlin);
  }
};

// FIXME remove ifdef when compiling in all the EOS
#ifdef PELEC_USE_SRK
template <>
struct CloseTransport<eos::SRK>
{
  void operator()()
  {
    amrex::The_Device_Arena()->free(trans_parm_g);

    amrex::The_Arena()->free(trans_parm.trans_wt);
    amrex::The_Arena()->free(trans_parm.trans_iwt);
    amrex::The_Arena()->free(trans_parm.trans_eps);
    amrex::The_Arena()->free(trans_parm.trans_sig);
    amrex::The_Arena()->free(trans_parm.trans_dip);
    amrex::The_Arena()->free(trans_parm.trans_pol);
    amrex::The_Arena()->free(trans_parm.trans_zrot);
    amrex::The_Arena()->free(trans_parm.trans_fitmu);
    amrex::The_Arena()->free(trans_parm.trans_fitlam);
    amrex::The_Arena()->free(trans_parm.trans_fitdbin);
    amrex::The_Arena()->free(trans_parm.trans_nlin);

    amrex::The_Arena()->free(trans_parm.Afac);
    amrex::The_Arena()->free(trans_parm.Bfac);
    amrex::The_Arena()->free(trans_parm.Kappai);
    amrex::The_Arena()->free(trans_parm.Sigmaij);
    amrex::The_Arena()->free(trans_parm.omega);
  }
};
#endif

template <typename EOSType>
struct SRKChungCorrections
{
  template <class... Args>
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(Args...)
  {
  }
};

// FIXME remove ifdef when compiling in all the EOS
#ifdef PELEC_USE_SRK
template <>
struct SRKChungCorrections<eos::SRK>
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    bool wtr_get_mu,
    bool wtr_get_lam,
    amrex::Real& Temp,
    amrex::Real* Xloc,
    amrex::Real& rholoc,
    amrex::Real& wbar,
    amrex::Real& mu,
    amrex::Real& xi,
    amrex::Real& lam,
    TransParm const* trans_parm)
  {
    amrex::Real sigma_M_3 = 0.0;
    amrex::Real Epsilon_M = 0.0;
    amrex::Real Omega_M = 0.0;
    amrex::Real MW_m = 0.0;
    amrex::Real DP_m_4 = 0.0;
    amrex::Real KappaM = 0.0;

    // Note: all the square roots could be precalculated

    for (int i = 0; i < NUM_SPECIES; ++i) {
      for (int j = 0; j < NUM_SPECIES; ++j) {
        amrex::Real T2 = trans_parm->trans_sig[i] * trans_parm->trans_sig[j];
        amrex::Real T3 = T2 * std::sqrt(T2);
        sigma_M_3 += Xloc[i] * Xloc[j] * T3;

        amrex::Real Epsilon_ij =
          std::sqrt(trans_parm->trans_eps[i] * trans_parm->trans_eps[j]);
        Epsilon_M += Xloc[i] * Xloc[j] * Epsilon_ij * T3;

        amrex::Real Omega_ij =
          0.5 * (trans_parm->omega[i] + trans_parm->omega[j]);
        Omega_M += Xloc[i] * Xloc[j] * Omega_ij * T3;

        amrex::Real MW_ij =
          2.0 / (trans_parm->trans_iwt[i] + trans_parm->trans_iwt[j]);
        MW_m += Xloc[i] * Xloc[j] * Epsilon_ij * T2 * std::sqrt(MW_ij);

        DP_m_4 += Xloc[i] * Xloc[j] * trans_parm->trans_dip[i] *
                  trans_parm->trans_dip[i] * trans_parm->trans_dip[j] *
                  trans_parm->trans_dip[j] / (T3 * Epsilon_ij);

        KappaM += Xloc[i] * Xloc[j] *
                  std::sqrt(trans_parm->Kappai[i] * trans_parm->Kappai[j]);
      }
    }

    MW_m *= MW_m;
    amrex::Real sigma_M = std::cbrt(sigma_M_3);
    amrex::Real InvSigma1 = 1.0 / sigma_M;
    amrex::Real InvSigma3 = 1.0 / sigma_M_3;
    Epsilon_M *= InvSigma3;
    amrex::Real Tstar = Temp / Epsilon_M;
    amrex::Real Tcm = 1.2593 * Epsilon_M;
    amrex::Real Vcm = 1.8887 * sigma_M_3;
    Omega_M *= InvSigma3;
    MW_m = MW_m * InvSigma3 * InvSigma1 / (Epsilon_M * Epsilon_M);
    DP_m_4 = DP_m_4 * sigma_M_3 * Epsilon_M;
    amrex::Real DP_red_4 = 297.2069113e6 * DP_m_4 / (Vcm * Vcm * Tcm * Tcm);
    amrex::Real y = Vcm * rholoc / (6.0 * wbar);
    amrex::Real A[10];
    for (int i = 0; i < 10; ++i) {
      A[i] = trans_parm->Afac[4 * i] + trans_parm->Afac[4 * i + 1] * Omega_M +
             trans_parm->Afac[4 * i + 2] * DP_red_4 +
             trans_parm->Afac[4 * i + 3] * KappaM;
    }
    amrex::Real G1 = (1.0 - 0.5 * y) / ((1.0 - y) * (1.0 - y) * (1.0 - y));
    amrex::Real G2 = (A[0] * (1.0 - std::exp(-A[3] * y)) / y +
                      A[1] * G1 * std::exp(A[4] * y) + A[2] * G1) /
                     (A[0] * A[3] + A[1] + A[2]);
    amrex::Real eta_P =
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
    amrex::Real beta = 1.0 / H2 + B[5] * y;

    // Set nonideal conductivity
    if (wtr_get_lam) {
      lam = lam * beta + lambda_p;
      // lam = H2;
      // mu = H2;
      // xi = beta;
    }
  }
};
#endif

template <typename EOSType>
struct BinaryDiff
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const amrex::Real wbar,
    amrex::Real* Xloc,
    amrex::Real* Yloc,
    amrex::Real* logT,
    amrex::Real& rholoc,
    amrex::Real& Tloc,
    amrex::Real* Ddiag,
    TransParm const* trans_parm)
  {
    amrex::Real term1, term2, dbintemp;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      term1 = 0.0;
      term2 = 0.0;
      for (int j = 0; j < NUM_SPECIES; ++j) {
        if (i != j) {
          dbintemp =
            trans_parm->trans_fitdbin[4 * (i + NUM_SPECIES * j)] +
            trans_parm->trans_fitdbin[1 + 4 * (i + NUM_SPECIES * j)] * logT[0] +
            trans_parm->trans_fitdbin[2 + 4 * (i + NUM_SPECIES * j)] * logT[1] +
            trans_parm->trans_fitdbin[3 + 4 * (i + NUM_SPECIES * j)] * logT[2];
          term1 = term1 + Yloc[j];
          term2 = term2 + Xloc[j] / std::exp(dbintemp);
        }
      }
      Ddiag[i] = trans_parm->trans_wt[i] * term1 / term2 / wbar;
    }

    // Call CKRP ?
    const amrex::Real pscale =
      Constants::PATM * wbar / (Constants::RU * Tloc * rholoc);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      Ddiag[i] = rholoc * pscale * Ddiag[i];
    }
  }
};

// FIXME remove ifdef when compiling in all the EOS
#ifdef PELEC_USE_SRK
template <>
struct BinaryDiff<eos::SRK>
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const amrex::Real wbar,
    amrex::Real* Xloc,
    amrex::Real* Yloc,
    amrex::Real* logT,
    amrex::Real& rholoc,
    amrex::Real& Tloc,
    amrex::Real* Ddiag,
    TransParm const* trans_parm)
  {
    amrex::Real term1, term2, Upsilonij;
    amrex::Real dbintemp;
    const amrex::Real Pst = 1013250.00;
    int idx_ij, idx_ik, idx_jk;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      term1 = 0.0;
      term2 = 0.0;
      for (int j = 0; j < NUM_SPECIES; ++j) {
        idx_ij = i + NUM_SPECIES * j;
        if (i != j) {
          dbintemp = trans_parm->trans_fitdbin[4 * (idx_ij)] +
                     trans_parm->trans_fitdbin[1 + 4 * (idx_ij)] * logT[0] +
                     trans_parm->trans_fitdbin[2 + 4 * (idx_ij)] * logT[1] +
                     trans_parm->trans_fitdbin[3 + 4 * (idx_ij)] * logT[2];
          dbintemp = std::exp(dbintemp);

          Upsilonij = 0.0;
          for (int k = 0; k < NUM_SPECIES; ++k) {
            idx_ik = i + NUM_SPECIES * k;
            idx_jk = j + NUM_SPECIES * k;
            Upsilonij +=
              trans_parm->trans_iwt[k] * Yloc[k] *
              (8.0 *
                 (trans_parm->Sigmaij[idx_ik] * trans_parm->Sigmaij[idx_ik] *
                    trans_parm->Sigmaij[idx_ik] +
                  trans_parm->Sigmaij[idx_jk] * trans_parm->Sigmaij[idx_jk] *
                    trans_parm->Sigmaij[idx_jk]) -
               6.0 *
                 (trans_parm->Sigmaij[idx_ik] * trans_parm->Sigmaij[idx_ik] +
                  trans_parm->Sigmaij[idx_jk] * trans_parm->Sigmaij[idx_jk]) *
                 trans_parm->Sigmaij[idx_ij] -
               3.0 *
                 std::pow(
                   (trans_parm->Sigmaij[idx_ik] * trans_parm->Sigmaij[idx_ik] -
                    trans_parm->Sigmaij[idx_jk] * trans_parm->Sigmaij[idx_jk]),
                   2.0) /
                 trans_parm->Sigmaij[idx_ij] +
               trans_parm->Sigmaij[idx_ij] * trans_parm->Sigmaij[idx_ij] *
                 trans_parm->Sigmaij[idx_ij]);
          }
          Upsilonij =
            Upsilonij * rholoc * Constants::Avna * 3.141592653589793 / 12.0 +
            1.0;
          dbintemp = dbintemp * Pst /
                     (Constants::RU * Tloc * Upsilonij); // * wbar/rholoc
          term1 += Yloc[j];
          term2 += Xloc[j] / dbintemp;
        }
      }
      Ddiag[i] = trans_parm->trans_wt[i] * term1 / term2; // * rholoc / wbar
    }
  }
};
#endif

struct SimpleTransport
{
  using transport_type = SimpleTransport;

  static std::string identifier() { return "SimpleTransport"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void comp_pure_bulk(
    amrex::Real Tloc,
    amrex::Real* muloc,
    amrex::Real* xiloc,
    TransParm const* trans_parm)
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
      if (trans_parm->trans_nlin[i] == 0) {
        cvkint[i] = 0.0;
        cvkrot[i] = 0.0;
      } else if (trans_parm->trans_nlin[i] == 1) {
        cvkint[i] = cvk[i] * trans_parm->trans_wt[i] / Constants::RU - 1.50;
        cvkrot[i] = 1.0;
      } else {
        cvkint[i] = cvk[i] * trans_parm->trans_wt[i] / Constants::RU - 1.50;
        cvkrot[i] = 1.5;
      }
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      epskoverTstd = trans_parm->trans_eps[i] / 298.0;
      epskoverT = trans_parm->trans_eps[i] / Tloc;
      Fnorm[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverTstd) +
                 (2.0 + .50 * M_PI * M_PI) * epskoverTstd +
                 pi3_2 * std::sqrt(epskoverTstd) * epskoverTstd;

      FofT[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverT) +
                (2.0 + .50 * M_PI * M_PI) * epskoverT +
                pi3_2 * std::sqrt(epskoverT) * epskoverT;
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      if (trans_parm->trans_nlin[i] == 0) {
        xiloc[i] = 0.0;
      } else {
        //   zrot/crot approximately zint / cint by assuming vibrational
        //   internal energy is small cvkrot is scaled by wk / Ru = mk / kb
        //   relative to standard specific cv
        xiloc[i] =
          0.250 * M_PI *
          (cvkint[i] / (cvkint[i] + 1.50) * cvkint[i] / (cvkint[i] + 1.50)) *
          trans_parm->trans_zrot[i] / cvkrot[i] * Fnorm[i] / FofT[i] * muloc[i];
      }
    }
  }

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
    amrex::Real trace = 1.e-15;
    // Call CKRP ?
    amrex::Real wbar;
    amrex::Real Xloc[NUM_SPECIES];
    amrex::Real muloc[NUM_SPECIES];
    amrex::Real xiloc[NUM_SPECIES];
    amrex::Real lamloc[NUM_SPECIES];
    amrex::Real logT[NUM_FIT - 1];

    logT[0] = std::log(Tloc);
    logT[1] = logT[0] * logT[0];
    logT[2] = logT[0] * logT[1];

    int nspec = NUM_SPECIES;

    amrex::Real sum = 0.;

    for (int i = 0; i < nspec; ++i) {
      sum = sum + Yloc[i];
    }

    wbar = 0.;
    amrex::Real real_nspec = NUM_SPECIES;

    for (int i = 0; i < NUM_SPECIES; ++i) {
      Yloc[i] = Yloc[i] + trace * (sum / real_nspec - Yloc[i]);
    }

    for (int i = 0; i < NUM_SPECIES; ++i) {
      wbar = wbar + Yloc[i] * trans_parm->trans_iwt[i];
    }
    wbar = 1.0 / wbar;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      Xloc[i] = Yloc[i] * wbar * trans_parm->trans_iwt[i];
    }
    if (wtr_get_mu) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        muloc[i] = trans_parm->trans_fitmu[4 * i] +
                   trans_parm->trans_fitmu[1 + 4 * i] * logT[0] +
                   trans_parm->trans_fitmu[2 + 4 * i] * logT[1] +
                   trans_parm->trans_fitmu[3 + 4 * i] * logT[2];
        muloc[i] = std::exp(muloc[i]);
      }

      mu = 0.0;

      for (int i = 0; i < NUM_SPECIES; ++i) {
        mu = mu + Xloc[i] * (muloc[i] * muloc[i] * muloc[i] * muloc[i] *
                             muloc[i] * muloc[i]);
      }
      mu = std::cbrt(std::sqrt(mu)); // mu = (sum[ Xloc_i * muloc_i^6])^(1/6)

      //  assumption that we only get bulk viscosity if we are already getting
      //  shear viscosity

      if (wtr_get_xi) {
        comp_pure_bulk(Tloc, muloc, xiloc, trans_parm);
        xi = 0.0;
        for (int i = 0; i < NUM_SPECIES; ++i) {
          xi = xi + Xloc[i] * (std::sqrt(std::sqrt(
                                xiloc[i] * xiloc[i] *
                                xiloc[i]))); // xi = sum[Xloc_i*xiloc_i^(3/4)];
        }
        xi = std::cbrt(
          xi * xi * xi * xi); // xi = (sum[Xloc_i*xiloc_i^(3/4)])^(4/3)
      }
    }

    if (wtr_get_lam) {
      for (int i = 0; i < NUM_SPECIES; ++i) {
        lamloc[i] = trans_parm->trans_fitlam[4 * i] +
                    trans_parm->trans_fitlam[1 + 4 * i] * logT[0] +
                    trans_parm->trans_fitlam[2 + 4 * i] * logT[1] +
                    trans_parm->trans_fitlam[3 + 4 * i] * logT[2];
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
    SRKChungCorrections<EosType>()(
      wtr_get_mu, wtr_get_lam, Tloc, Xloc, rholoc, wbar, mu, xi, lam,
      trans_parm);

    if (wtr_get_Ddiag) {
      BinaryDiff<EosType>()(
        wbar, Xloc, Yloc, logT, rholoc, Tloc, Ddiag, trans_parm);
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
  AMREX_GPU_HOST_DEVICE SimpleTransport(Args...)
  {
  }
};
} // namespace transport
using TransportType = transport::SimpleTransport;
} // namespace physics
} // namespace pele
#endif
