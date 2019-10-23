#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

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
