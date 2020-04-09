#include "chemistry_file.H"

void CKSYMS(int * kname, int * plenkname) {}

/* Returns R, Rc, Patm */
void CKRP(double *  ru, double *  ruc, double *  pa)
{
     *ru  = 8.31451e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}

AMREX_GPU_HOST_DEVICE void CKINIT()
{
}

AMREX_GPU_HOST_DEVICE void CKFINALIZE()
{
}

