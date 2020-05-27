#ifndef _main_F_H_
#define _main_F_H_

#include <AMReX_BLFort.H>

#include <AMReX.H>

extern "C"
{
    void initialize_data
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(temperature),
     const BL_FORT_FAB_ARG_3D(density),
     const BL_FORT_FAB_ARG_3D(energy),
     const amrex::Real* dx, const amrex::Real* plo, const amrex::Real* phi);

}

#endif // _main_F_H_
