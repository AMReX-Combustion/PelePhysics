#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
initialize_data(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& energy,
		amrex::GeometryData const& geomdata) noexcept
{
   const amrex::Real* plo = geomdata.ProbLo();
   const amrex::Real* phi = geomdata.ProbHi();
   const amrex::Real* dx  = geomdata.CellSize();

	amrex::Real dTemp = 5.0;
	amrex::Real dRho = 0.005;
	amrex::Real denergy = 100.0;
	amrex::Real y = plo[1] + (j+0.5)*dx[1];
	amrex::Real x = plo[0] + (i+0.5)*dx[0];
	amrex::Real pi = 3.1415926535897932;
	amrex::GpuArray<amrex::Real,3> L;
	amrex::GpuArray<amrex::Real,3> P;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> Y_lo;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> Y_hi;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> ei;

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
#if ( AMREX_SPACEDIM == 1 )
	temp(i,j,k) = 500.0;
	rho(i,j,k) = 0.01;
#else
	temp(i,j,k) = 500.0 + dTemp * std::sin(2.0*pi*y/P[1]);
	rho(i,j,k) = 0.01 + dRho * std::sin(2.0*pi*y/P[1]);
#endif
	for (int n = 0; n < NUM_SPECIES; n++) {
            mf(i,j,k,n) = Y_lo[n] + (Y_hi[n]-Y_lo[n]) * x / L[0];
	}
        // corr Yk
	amrex::Real dummy = 0.0;
	for (int n = 0; n < NUM_SPECIES-1; n++) {
	    dummy = dummy + mf(i,j,k,n);
	}
	mf(i,j,k,NUM_SPECIES-1) = 1.0 - dummy;
        // compute real E + perturb it for future tests
        auto eos = pele::physics::PhysicsType::eos();
        eos.T2Ei(temp(i,j,k), ei.arr);
	energy(i,j,k) = 0.0;
	for (int n = 0; n < NUM_SPECIES; n++) {
	    energy(i,j,k) = energy(i,j,k) + mf(i,j,k,n) * ei[n];
	}
#if ( AMREX_SPACEDIM == 1 )
	energy(i,j,k) = energy(i,j,k); 
#else
	energy(i,j,k) = energy(i,j,k) + denergy * std::sin(2.0*pi*y/P[1]); 
#endif
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cp(int i, int j, int k,  
	amrex::Array4<const amrex::Real> const& mf,
	amrex::Array4<const amrex::Real> const& temp,
	amrex::Array4<amrex::Real> const& cp) noexcept 
{

	amrex::GpuArray<amrex::Real,NUM_SPECIES> mf_pt;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}

        auto eos = pele::physics::PhysicsType::eos();
        eos.TY2Cp(temp(i,j,k), mf_pt.arr, cp(i,j,k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cv(int i, int j, int k,  
	amrex::Array4<const amrex::Real> const& mf,
	amrex::Array4<const amrex::Real> const& temp,
	amrex::Array4<amrex::Real> const& cv) noexcept 
{

	amrex::GpuArray<amrex::Real,NUM_SPECIES> mf_pt;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}

        auto eos = pele::physics::PhysicsType::eos();
	eos.TY2Cv(temp(i,j,k), mf_pt.arr, cv(i,j,k));
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_T_from_EY(int i, int j, int k,  
	amrex::Array4<const amrex::Real> const& mf,
	amrex::Array4<amrex::Real> const& temp,
	amrex::Array4<const amrex::Real> const& energy) noexcept 
{

	amrex::GpuArray<amrex::Real,NUM_SPECIES> mf_pt;
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}

        auto eos = pele::physics::PhysicsType::eos();
	eos.EY2T(energy(i,j,k), mf_pt.arr, temp(i,j,k));
}


#endif
