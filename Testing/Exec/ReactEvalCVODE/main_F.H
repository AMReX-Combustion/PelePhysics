#ifndef _main_F_H_
#define _main_F_H_

#include <AMReX_BLFort.H>

extern "C"
{
    void extern_init
    (const int* name, const int* namlen, const int* cvode_iE);

    void extern_close();

    void get_num_spec(int* nspecies);

    void initialize_data
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(rY_src_ext),
     const BL_FORT_FAB_ARG_3D(energy),
     const BL_FORT_FAB_ARG_3D(energy_src_ext),
     const Real* dx, const Real* plo, const Real* phi);

}

#endif // _main_F_H_
