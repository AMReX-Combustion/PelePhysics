#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include "mechanism.h"

#include <EOS.H>

AMREX_GPU_HOST_DEVICE
inline
void
initialize_data(int i, int j, int k, int fuel_id, 
        amrex::Array4<amrex::Real> const& rhoY,
        amrex::Array4<amrex::Real> const& frcExt,
        amrex::Array4<amrex::Real> const& rhoE,
        amrex::Array4<amrex::Real> const& frcEExt,
        GpuArray<Real, AMREX_SPACEDIM>  dx, 
        GpuArray<Real, AMREX_SPACEDIM>  plo, 
        GpuArray<Real, AMREX_SPACEDIM>  phi ) noexcept
{
    amrex::Real Temp_lo = 2000.0;
    amrex::Real Temp_hi = 2500.0;
    amrex::Real dTemp = 100.0;
    amrex::Real pressure = 1013250.0;
    amrex::Real density, energy, temp;
    amrex::GpuArray<amrex::Real,NUM_SPECIES> X;
    amrex::GpuArray<amrex::Real,NUM_SPECIES> Y;
    amrex::Real y = plo[1] + (j+0.5)*dx[1];
    amrex::Real x = plo[0] + (i+0.5)*dx[0];
    amrex::Real pi = 3.1415926535897932;
    amrex::GpuArray<amrex::Real,3> L;
    amrex::GpuArray<amrex::Real,3> P;

    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        L[n] = phi[n] - plo[n];
        P[n] = L[n] / 4.0;
    }
        // Y
    for (int n = 0; n < NUM_SPECIES; n++) {
        X[n] = 0.0;
    }
    X[O2_ID]   = 0.2;
    X[fuel_id] = 0.1;
    X[N2_ID]   = 0.7;
    EOS::X2Y(&X[0],&Y[0]);
    // T
    temp =  Temp_lo + (Temp_hi-Temp_lo)*y/L[1] + dTemp * std::sin(2.0*pi*y/P[1]);
    // get rho and E 
    EOS::PYT2RE(pressure, &Y[0], temp, density, energy);
    // Fill vect
    for (int n = 0; n < NUM_SPECIES; n++) {
        rhoY(i,j,k,n) = Y[n]*density;   
        frcExt(i,j,k,n) = 0.0;
    }
    rhoY(i,j,k,NUM_SPECIES) = temp; 
#if defined(USE_SUNDIALS_PP) || defined(USE_RK64_PP)
    rhoE(i,j,k) = energy * density; 
#else
    rhoE(i,j,k) = energy; 
#endif
    frcEExt(i,j,k) = 0.0;
}


AMREX_GPU_HOST_DEVICE
inline
void
gpu_flatten(int icell, int i, int j, int k,  
        amrex::Array4<amrex::Real> const& rhoY,
        amrex::Array4<amrex::Real> const& frcExt,
        amrex::Array4<amrex::Real> const& rhoE,
        amrex::Array4<amrex::Real> const& frcEExt,
        amrex::Real * tmp_vect, amrex::Real * tmp_src_vect,
        amrex::Real * tmp_vect_energy, amrex::Real * tmp_src_vect_energy) noexcept
{
    amrex::Real * TV_cell  = tmp_vect + icell * (NUM_SPECIES+1);
    amrex::Real * TVS_cell = tmp_src_vect + icell * (NUM_SPECIES);
    for (int n = 0; n < NUM_SPECIES; n++) {
        TV_cell[n]   = rhoY(i,j,k,n);
        TVS_cell[n]  = frcExt(i,j,k,n);
    }
    TV_cell[NUM_SPECIES]       = rhoY(i,j,k,NUM_SPECIES);  
    tmp_vect_energy[icell]     = rhoE(i,j,k,0);
    tmp_src_vect_energy[icell] = frcEExt(i,j,k,0); 
}


AMREX_GPU_HOST_DEVICE
inline
void
gpu_unflatten(int icell, int i, int j, int k,  
        amrex::Array4<amrex::Real> const& rhoY,
        amrex::Array4<amrex::Real> const& rhoE,
        amrex::Real * tmp_vect,
        amrex::Real * tmp_vect_energy) noexcept
{
    amrex::Real * TV_cell  = tmp_vect + icell * (NUM_SPECIES+1);
    for (int n = 0; n < NUM_SPECIES; n++) {
        rhoY(i,j,k,n) = TV_cell[n];
    }
    rhoY(i,j,k,NUM_SPECIES) = TV_cell[NUM_SPECIES];  
    rhoE(i,j,k,0) = tmp_vect_energy[icell];
}




#endif
