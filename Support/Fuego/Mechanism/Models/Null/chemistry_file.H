#include <AMReX_Gpu.H>

extern "C"
{
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa);
AMREX_GPU_HOST_DEVICE void CKINIT();
AMREX_GPU_HOST_DEVICE void CKFINALIZE();
AMREX_GPU_HOST_DEVICE void CKWC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  wdot);
}
