#ifndef _TRANSPORT_PARAMS_H_
#define _TRANSPORT_PARAMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace transport_params {

void init();
void finalize();

extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_wt;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_iwt;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_eps;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_sig;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_dip;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_pol;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_zrot;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_fitmu;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_fitlam;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* trans_fitdbin;
extern AMREX_GPU_DEVICE_MANAGED int* trans_nlin;
extern AMREX_GPU_DEVICE_MANAGED int array_size;
extern AMREX_GPU_DEVICE_MANAGED int fit_length;

} // namespace transport_params

#endif
