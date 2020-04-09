#ifndef _CHEMISTRY_FILE_H_
#define _CHEMISTRY_FILE_H_
void CKSYMS(int * kname, int * lenkname);
void CKRP(double *  ru, double *  ruc, double *  pa);
AMREX_GPU_HOST_DEVICE void CKINIT();
AMREX_GPU_HOST_DEVICE void CKFINALIZE();
#endif
