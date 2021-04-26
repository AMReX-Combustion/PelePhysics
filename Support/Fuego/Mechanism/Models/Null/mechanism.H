#ifndef MECHANISM_H
#define MECHANISM_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  X
*/
#endif

/* Species */
#define XO_ID 0

#define NUM_ELEMENTS 1
#define NUM_SPECIES 1
#define NUM_REACTIONS 0

/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
/* INIT and FINALIZE  */
void CKINIT();
void CKFINALIZE();
/* MISC  */
void CKSYMS_STR(amrex::Vector<std::string>& kname);


/* Returns R, Rc, Patm */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}

/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
     wdot[0] = 0.0e0;
}

#endif

/* End of file  */
