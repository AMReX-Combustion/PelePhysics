#ifndef _kernel_H_
#define _kernel_H_
#include <AMReX_GpuDevice.H>

#define NUM_SPECIES 21
#define NUM_REACTIONS 84

using namespace amrex;

AMREX_GPU_DEVICE void W_spec_d(Real rho,
                               Real temp,
                               Real * Y, int strideY,
                               Real * wdot);
void kinit();

#endif
