#include <cstdio>

#include <AMReX_Arena.H>

#include "mechanism.h"
#include "chemistry_file.H"
#include "TransportParams.H"

namespace transport_params {

AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_wt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_iwt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_eps;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_sig;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_dip;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_pol;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_zrot;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_fitmu;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_fitlam;
AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_fitdbin;
AMREX_GPU_DEVICE_MANAGED int* trans_nlin;
AMREX_GPU_DEVICE_MANAGED int array_size = NUM_SPECIES;
AMREX_GPU_DEVICE_MANAGED int fit_length = NUM_FIT;

void
init()
{
  //    std::cout << " array_size " << array_size << std::endl;
  //    std::cout << " fit_length " << fit_length << std::endl;
  trans_wt = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  trans_iwt = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  trans_eps = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  trans_sig = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  dip = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  trans_pol = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  trans_zrot = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));

  trans_fitmu = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * fit_length));
  trans_fitlam = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * fit_length));
  trans_fitdbin = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * array_size * fit_length));

  trans_nlin = static_cast<int*>(
    amrex::The_Managed_Arena()->alloc(sizeof(int) * array_size));

  egtransetWT(trans_wt);
  egtransetEPS(trans_eps);
  egtransetSIG(trans_sig);
  egtransetDIP(trans_dip);
  egtransetPOL(trans_pol);
  egtransetZROT(trans_zrot);
  egtransetNLIN(trans_nlin);
  egtransetCOFETA(trans_fitmu);
  egtransetCOFLAM(trans_fitlam);
  egtransetCOFD(trans_fitdbin);

  for (int i = 0; i < array_size; ++i) {
    trans_iwt[i] = 1. / trans_wt[i];
  }
}

void
finalize()
{
  amrex::The_Managed_Arena()->free(trans_wt);
  amrex::The_Managed_Arena()->free(trans_iwt);
  amrex::The_Managed_Arena()->free(trans_eps);
  amrex::The_Managed_Arena()->free(trans_sig);
  amrex::The_Managed_Arena()->free(trans_dip);
  amrex::The_Managed_Arena()->free(trans_pol);
  amrex::The_Managed_Arena()->free(trans_zrot);
  amrex::The_Managed_Arena()->free(trans_fitmu);
  amrex::The_Managed_Arena()->free(trans_fitlam);
  amrex::The_Managed_Arena()->free(trans_fitdbin);
  amrex::The_Managed_Arena()->free(trans_nlin);
}

} // namespace transport_params
