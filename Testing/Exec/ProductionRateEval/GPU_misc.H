#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <EOS.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
inline
void
initialize_data(int i, int j, int k,  
		Array4<Real> const& mf,
		Array4<Real> const& temp,
		Array4<Real> const& rho,
                std::vector<Real> dx, 
		std::vector<Real> plo, 
		std::vector<Real> phi ) noexcept
{
	Real dTemp = 5.0;
	Real dRho = 0.005;
	Real denergy = 100.0;
	Real y = plo[1] + (j+0.5)*dx[1];
	Real x = plo[0] + (i+0.5)*dx[0];
	Real pi = 3.1415926535897932;
	GpuArray<Real,3> L;
	GpuArray<Real,3> P;
	GpuArray<Real,NUM_SPECIES> Y_lo;
	GpuArray<Real,NUM_SPECIES> Y_hi;

	for (int n = 0; n < 3; n++) {
		L[n] = phi[n] - plo[n];
		P[n] = L[n] / 4.0;
	}
	for (int n = 0; n < NUM_SPECIES; n++) {
		Y_lo[n] = 0.0;
		Y_hi[n] = 1.0 / NUM_SPECIES ;
	}
	Y_lo[0] = 1.0;

	// T, Yk, rho
        Real Tavg = 500;
        Real Ravg = 0.01;
#if ( AMREX_SPACEDIM == 1 )
	temp(i,j,k) = Tavg;
	rho(i,j,k) = Ravg;
#else
	temp(i,j,k) = Tavg + dTemp * std::sin(2.0*pi*y/P[1]);
	rho(i,j,k) = Ravg + dRho * std::sin(2.0*pi*y/P[1]);
#endif
	for (int n = 0; n < NUM_SPECIES; n++) {
            mf(i,j,k,n) = Y_lo[n] + (Y_hi[n]-Y_lo[n]) * x / L[0];
	}
        // corr Yk
	Real dummy = 0.0;
	for (int n = 0; n < NUM_SPECIES-1; n++) {
	    dummy = dummy + mf(i,j,k,n);
	}
	mf(i,j,k,NUM_SPECIES-1) = 1.0 - dummy;
}


AMREX_GPU_HOST_DEVICE
inline
void
get_cp(int i, int j, int k,  
	Array4<Real> const& mf,
	Array4<Real> const& temp,
	Array4<Real> const& cp) noexcept 
{

	GpuArray<Real,NUM_SPECIES> mf_pt;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}

	EOS::TY2Cp(temp(i,j,k), mf_pt.arr, cp(i,j,k));
}

AMREX_GPU_HOST_DEVICE
inline
void
get_cv(int i, int j, int k,  
	Array4<Real> const& mf,
	Array4<Real> const& temp,
	Array4<Real> const& cv) noexcept 
{

	GpuArray<Real,NUM_SPECIES> mf_pt;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}

	EOS::TY2Cv(temp(i,j,k), mf_pt.arr, cv(i,j,k));
}


AMREX_GPU_HOST_DEVICE
inline
void
get_T_from_EY(int i, int j, int k,  
	Array4<Real> const& mf,
	Array4<Real> const& temp,
	Array4<Real> const& energy) noexcept 
{

	GpuArray<Real,NUM_SPECIES> mf_pt;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}

	EOS::EY2T(energy(i,j,k), mf_pt.arr, temp(i,j,k));
}


#endif
