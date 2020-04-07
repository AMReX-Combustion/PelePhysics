#ifndef _Transport_F_H_
#define _Transport_F_H_
#include <AMReX_BLFort.H>

#include <AMReX.H>

extern "C"
{
  void get_transport_coeffs
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(temperature),
     const BL_FORT_FAB_ARG_3D(density),
     const BL_FORT_FAB_ARG_3D(rhoD),
     const BL_FORT_FAB_ARG_3D(mu),
     const BL_FORT_FAB_ARG_3D(xi),
     const BL_FORT_FAB_ARG_3D(lambda));
  
  void get_transport_coeffs_aux
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(massfraction),
     const BL_FORT_FAB_ARG_3D(temperature),
     const BL_FORT_FAB_ARG_3D(density),
     const BL_FORT_FAB_ARG_3D(rhoDaux));
}

#endif // _Transport_F_H_
