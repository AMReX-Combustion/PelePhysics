#ifndef _kernel_H_
#define _kernel_H_
#include <AMReX_GpuDevice.H>

#define NUM_SPECIES 9
#define NUM_REACTIONS 84

using namespace amrex;

AMREX_GPU_DEVICE Real Q_reac_d(Real rho,
                               Real temp,
                               Real * Y, int strideY,
                               int reacIdx);
void kinit();

#endif
