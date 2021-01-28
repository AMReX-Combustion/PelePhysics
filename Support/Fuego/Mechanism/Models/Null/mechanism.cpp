#include "chemistry_file.H"

void CKSYMS_STR(amrex::Vector<std::string>& kname) {}

/* Returns R, Rc, Patm */
void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}

AMREX_GPU_HOST_DEVICE void CKINIT()
{
}

AMREX_GPU_HOST_DEVICE void CKFINALIZE()
{
}

AMREX_GPU_HOST_DEVICE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
     wdot[0] = 0.0e0;
}

