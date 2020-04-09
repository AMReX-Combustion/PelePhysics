#ifndef _TRANSPORT_H_
#define _TRANSPORT_H_

//#include <limits>
//#include <cmath>

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>

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
comp_pure_bulk(amrex::Real Tloc, amrex::Real* muloc, amrex::Real* xiloc)
{
  amrex::Real cvk[NUM_SPECIES];
  amrex::Real cvkint[NUM_SPECIES];
  amrex::Real cvkrot[NUM_SPECIES];
  amrex::Real FofT[NUM_SPECIES];
  amrex::Real Fnorm[NUM_SPECIES];
  amrex::Real epskoverT, epskoverTstd;
  amrex::Real pi = 3.141592653589793238;
  amrex::Real pi3_2 = 5.56832799683171;
  amrex::Real Ru = 8.31446261815324e7;;

  CKCVMS(&Tloc, cvk);

  for (int i = 0; i < NUM_SPECIES; ++i) {
    if (trans_nlin[i] == 0) {
      cvkint[i] = 0.0;
      cvkrot[i] = 0.0;
    } else if (trans_nlin[i] == 1) {
      cvkint[i] = cvk[i] * trans_wt[i] / Ru - 1.50;
      cvkrot[i] = 1.0;
    } else {
      cvkint[i] = cvk[i] * wt[i] / Ru - 1.50;
      cvkrot[i] = 1.5;
    }
  }

  for (int i = 0; i < NUM_SPECIES; ++i) {
    epskoverTstd = trans_eps[i] / 298.0;
    epskoverT = trans_eps[i] / Tloc;
    Fnorm[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverTstd) +
               (2.0 + .50 * pi * pi) * epskoverTstd +
               pi3_2 * std::sqrt(epskoverTstd) * epskoverTstd;

    FofT[i] = 1.0 + 0.50 * pi3_2 * std::sqrt(epskoverT) +
              (2.0 + .50 * pi * pi) * epskoverT +
              pi3_2 * std::sqrt(epskoverT) * epskoverT;
  }

  for (int i = 0; i < NUM_SPECIES; ++i) {
    if (trans_nlin[i] == 0) {
      xiloc[i] = 0.0;
    } else {
      //   zrot/crot approximately zint / cint by assuming vibrational internal
      //   energy is small cvkrot is scaled by wk / Ru = mk / kb relative to
      //   standard specific cv
      xiloc[i] =
        0.250 * pi *
        (cvkint[i] / (cvkint[i] + 1.50) * cvkint[i] / (cvkint[i] + 1.50)) *
        trans_zrot[i] / cvkrot[i] * Fnorm[i] / FofT[i] * muloc[i];
    }
  }
}

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
  //  need to set Ru

  amrex::Real trace = 1.e-15;
  amrex::Real Patm = 1.01325e6;
  amrex::Real wbar, pscale;
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
    wbar = wbar + Yloc[i] * trans_iwt[i];
  }
  wbar = 1.0 / wbar;
  for (int i = 0; i < NUM_SPECIES; ++i) {
    Xloc[i] = Yloc[i] * wbar * trans_iwt[i];
  }
  if (wtr_get_mu) {
    for (int i = 0; i < NUM_SPECIES; ++i) {
      muloc[i] = trans_fitmu[4 * i] + trans_fitmu[1 + 4 * i] * logT[0] +
                 trans_fitmu[2 + 4 * i] * logT[1] + trans_fitmu[3 + 4 * i] * logT[2];
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
      comp_pure_bulk(Tloc, muloc, xiloc);
      xi = 0.0;
      for (int i = 0; i < NUM_SPECIES; ++i) {
        xi = xi + Xloc[i] * (std::sqrt(std::sqrt(
                              xiloc[i] * xiloc[i] *
                              xiloc[i]))); // xi = sum[Xloc_i*xiloc_i^(3/4)];
      }
      xi =
        std::cbrt(xi * xi * xi * xi); // xi = (sum[Xloc_i*xiloc_i^(3/4)])^(4/3)
    }
  }

  if (wtr_get_lam) {
    for (int i = 0; i < NUM_SPECIES; ++i) {
      lamloc[i] = trans_fitlam[4 * i] + trans_fitlam[1 + 4 * i] * logT[0] +
                  trans_fitlam[2 + 4 * i] * logT[1] + trans_fitlam[3 + 4 * i] * logT[2];
      lamloc[i] = std::exp(lamloc[i]);
    }

    lam = 0.;

    for (int i = 0; i < NUM_SPECIES; ++i) {
      lam +=
        Xloc[i] *
        std::sqrt(std::sqrt(lamloc[i])); // lam = sum[Xloc_i * (lamloc_i)^(1/4)]
    }
    lam = lam * lam * lam * lam; // lam = lam^4
  }

  if (wtr_get_Ddiag) {
    //       for (int i = 0; i < NUM_SPECIES ; ++i){
    //          for (int j = 0; j < i-1 ; ++j){
    //               dbinloc[i+NUM_SPECIES*j] =
    //               fitdbin[i+NUM_SPECIES*j]+fitdbin[1+4*(i+NUM_SPECIES*j)]*logT[1]
    //                   + fitdbin[2+4*(i+NUM_SPECIES*j)]*logT[2]+
    //                   fitdbin[3+4*(i+NUM_SPECIES*j)]*logT[3];
    //               dbinloc[i+NUM_SPECIES*j] =
    //               std::exp(dbinloc[i+NUM_SPECIES*j]);
    //           }
    //          dbinloc(i+NUM_SPECIES*i) = 0.0;
    //       }

    amrex::Real term1, term2, dbintemp;

    for (int i = 0; i < NUM_SPECIES; ++i) {
      term1 = 0.0;
      term2 = 0.0;
      for (int j = 0; j < NUM_SPECIES; ++j) {
        if (i != j) {
          dbintemp = trans_fitdbin[4 * (i + NUM_SPECIES * j)] +
                     trans_fitdbin[1 + 4 * (i + NUM_SPECIES * j)] * logT[0] +
                     trans_fitdbin[2 + 4 * (i + NUM_SPECIES * j)] * logT[1] +
                     trans_fitdbin[3 + 4 * (i + NUM_SPECIES * j)] * logT[2];
          term1 = term1 + Yloc[j];
          term2 = term2 + Xloc[j] / std::exp(dbintemp);
        }
      }
      Ddiag[i] = trans_wt[i] * term1 / term2 / wbar;
    }

    amrex::Real Ru = 8.31446261815324e7;;
    pscale = Patm * wbar / (Ru * Tloc * rholoc);

    for (int i = 0; i < NUM_SPECIES; ++i) {
      Ddiag[i] = rholoc * pscale * Ddiag[i];
    }
  }
}

#endif
