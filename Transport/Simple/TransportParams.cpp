#include <cstdio>

#include <AMReX_Arena.H>

#include "mechanism.h"
#include "chemistry_file.H"
#include "TransportParams.H"

namespace transport_params {

AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_wt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_iwt;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_eps;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_sig;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_dip;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_pol;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_zrot;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_fitmu;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_fitlam;
AMREX_GPU_DEVICE_MANAGED amrex::Real* eos_fitdbin;
AMREX_GPU_DEVICE_MANAGED int* eos_nlin;
AMREX_GPU_DEVICE_MANAGED int array_size;
AMREX_GPU_DEVICE_MANAGED int fit_length;

void
init()
{
  array_size = NUM_SPECIES;
  fit_length = NUM_FIT;
  //    std::cout << " array_size " << array_size << std::endl;
  //    std::cout << " fit_length " << fit_length << std::endl;
  eos_wt = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  eos_iwt = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  eos_eps = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  eos_sig = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  dip = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  eos_pol = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));
  eos_zrot = static_cast<amrex::Real*>(
    amrex::The_Managed_Arena()->alloc(sizeof(amrex::Real) * array_size));

  eos_fitmu = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * fit_length));
  eos_fitlam = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * fit_length));
  eos_fitdbin = static_cast<amrex::Real*>(amrex::The_Managed_Arena()->alloc(
    sizeof(amrex::Real) * array_size * array_size * fit_length));

  eos_nlin = static_cast<int*>(
    amrex::The_Managed_Arena()->alloc(sizeof(int) * array_size));

  egtransetWT(eos_wt);
  egtransetEPS(eos_eps);
  egtransetSIG(eos_sig);
  egtransetDIP(eos_dip);
  egtransetPOL(eos_pol);
  egtransetZROT(eos_zrot);
  egtransetNLIN(eos_nlin);
  egtransetCOFETA(eos_fitmu);
  egtransetCOFLAM(eos_fitlam);
  egtransetCOFD(eos_fitdbin);

  for (int i = 0; i < array_size; ++i) {
    eos_iwt[i] = 1. / eos_wt[i];
  }
}

void
finalize()
{
  amrex::The_Managed_Arena()->free(eos_wt);
  amrex::The_Managed_Arena()->free(eos_iwt);
  amrex::The_Managed_Arena()->free(eos_eps);
  amrex::The_Managed_Arena()->free(eos_sig);
  amrex::The_Managed_Arena()->free(eos_dip);
  amrex::The_Managed_Arena()->free(eos_pol);
  amrex::The_Managed_Arena()->free(eos_zrot);
  amrex::The_Managed_Arena()->free(eos_fitmu);
  amrex::The_Managed_Arena()->free(eos_fitlam);
  amrex::The_Managed_Arena()->free(eos_fitdbin);
  amrex::The_Managed_Arena()->free(eos_nlin);
}

} // namespace transport_params
