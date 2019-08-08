#ifndef _main_F_H_
#define _main_F_H_

#include <AMReX_BLFort.H>

#include <AMReX.H>

extern "C"
{
    void extern_init
    (const int* name, const int* namlen);

    void extern_close();
    
    void get_num_spec(int* nspecies);

    void initialize_data
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(temperature),
     const amrex::Real* dx, const amrex::Real* plo, const amrex::Real* phi);

    void get_cp
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(temperature),
     BL_FORT_FAB_ARG_3D(cp));
}

#endif // _main_F_H_
