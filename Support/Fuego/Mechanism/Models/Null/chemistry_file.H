#include <AMReX_Gpu.H>

extern "C"
{
void CKSYMS(int * kname, int * lenkname);
void CKRP(double *  ru, double *  ruc, double *  pa);
AMREX_GPU_HOST_DEVICE void CKINIT();
AMREX_GPU_HOST_DEVICE void CKFINALIZE();
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C, double *  wdot);
}
