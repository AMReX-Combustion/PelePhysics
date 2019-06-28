#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <Fuego_EOS.H>

AMREX_GPU_DEVICE
inline
void
gpu_RTY2W(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& cdots) noexcept
{
	amrex::GpuArray<amrex::Real,NUM_SPECIES> mf_pt;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> cdots_pts;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
		cdots_pts[n] = cdots(i,j,k,n);

	}

	EOS eos; 
	eos.eos_RTY2W(rho(i,j,k), temp(i,j,k), mf_pt.arr, cdots_pts.arr);
	//eos.eos_RTY2W(rho(i,j,k), temp(i,j,k), mf_pt.data(), cdots_pts.data());
	
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf(i,j,k,n) = mf_pt[n];
		cdots(i,j,k,n) = cdots_pts[n];

	}
}


#endif
