#ifndef AMREX_MISC_H_
#define AMREX_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include "mechanism.h"

#include <EOS.H>


AMREX_GPU_HOST_DEVICE
inline
void
box_flatten(int icell, int i, int j, int k, int ireactor_type,
            amrex::Array4<amrex::Real> const& rhoY,
            amrex::Array4<amrex::Real> const& frcExt,
            amrex::Array4<amrex::Real> const& temperature,
            amrex::Array4<amrex::Real> const& rhoE,
            amrex::Array4<amrex::Real> const& frcEExt,
            amrex::Array4<int> const& mask,
            amrex::Real* tmp_vect,
            amrex::Real* tmp_src_vect,
            amrex::Real* tmp_vect_energy,
            amrex::Real* tmp_src_vect_energy,
            int* tmp_mask_vect) noexcept
{
    amrex::Real mass_frac[NUM_SPECIES];
    amrex::Real rho = 0.0;
    amrex::Real Enrg_loc;
    amrex::Real temp;

    int offset_yvec = icell * (NUM_SPECIES+1);
    int offset_spec = icell * (NUM_SPECIES);

    for (int n = 0; n < NUM_SPECIES; n++) {
        tmp_vect[offset_yvec + n]     = rhoY(i,j,k,n);
        tmp_src_vect[offset_spec + n] = frcExt(i,j,k,n);
        rho += tmp_vect[offset_yvec + n];
    }
    amrex::Real rho_inv        = 1.0 / rho;
    temp                       = temperature(i,j,k);
    tmp_vect_energy[icell]     = rhoE(i,j,k);
    tmp_src_vect_energy[icell] = frcEExt(i,j,k);
    tmp_mask_vect[i]           = mask(i,j,k);

    /* T update with energy and Y */
    for (int n = 0; n < NUM_SPECIES; n++) {
        mass_frac[n] = tmp_vect[offset_yvec + n] * rho_inv;
    }
    Enrg_loc = tmp_vect_energy[icell] / rho;
    if (ireactor_type == 1){
        EOS::EY2T(Enrg_loc,mass_frac,temp);
    } else {
        EOS::HY2T(Enrg_loc,mass_frac,temp);
    }
    tmp_vect[offset_yvec + NUM_SPECIES] = temp;
}


AMREX_GPU_HOST_DEVICE
inline
void
box_unflatten(int icell, int i, int j, int k, int ireactor_type,
              amrex::Array4<amrex::Real> const& rhoY,
              amrex::Array4<amrex::Real> const& T_in,
              amrex::Array4<amrex::Real> const& rhoE,
              amrex::Array4<amrex::Real> const& frcEExt,
              amrex::Array4<amrex::Real> const& FC_in,
              amrex::Real* tmp_vect,
              amrex::Real* tmp_vect_energy,
              int* tmp_FCunt,
              amrex::Real dt) noexcept
{
    amrex::Real mass_frac[NUM_SPECIES];
    amrex::Real rho = 0.0;
    amrex::Real Enrg_loc;
    amrex::Real temp;

    int offset_yvec = icell * (NUM_SPECIES+1);

    for (int n = 0; n < NUM_SPECIES; n++) {
        rhoY(i,j,k,n) = tmp_vect[offset_yvec + n];
        rho += tmp_vect[offset_yvec + n];
    }
    amrex::Real rho_inv = 1.0 / rho;
    temp                = tmp_vect[offset_yvec + NUM_SPECIES];
    rhoE(i,j,k)         = tmp_vect_energy[icell] + dt * frcEExt(i,j,k,0);    

    /* T update with energy and Y */
    for (int n = 0; n < NUM_SPECIES; n++) {
        mass_frac[n] = rhoY(i,j,k,n) * rho_inv;
    }
    Enrg_loc = rhoE(i,j,k) / rho;
    if (ireactor_type == 1){
        EOS::EY2T(Enrg_loc,mass_frac,temp);
    } else {
        EOS::HY2T(Enrg_loc,mass_frac,temp);
    }
    T_in(i,j,k) = temp;

    FC_in(i,j,k) = tmp_FCunt[icell];
}
#endif
