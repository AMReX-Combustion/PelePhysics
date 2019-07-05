#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <Fuego_EOS.H>

AMREX_GPU_HOST_DEVICE
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

        //printf(" gpu_RTY2W \n");
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
		cdots_pts[n] = cdots(i,j,k,n);
	}
        
        //printf(" gpu_RTY2W %d %d %d  \n", i,j,k);
        //printf(" gpu_RTY2W %14.6e %14.6e \n", rho(i,j,k), temp(i,j,k));

	EOS eos; 
	eos.eos_RTY2W(rho(i,j,k), temp(i,j,k), mf_pt.arr, cdots_pts.arr);
	
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf(i,j,k,n) = mf_pt[n];
		cdots(i,j,k,n) = cdots_pts[n];

	}
}


__global__ void MyLaunchTest(int ncells,
                              amrex::Array4<amrex::Real> const rho,
                              amrex::Array4<amrex::Real> const temp,
                              amrex::Array4<amrex::Real> const mf,
                              amrex::Array4<amrex::Real> const cdots,
                              int lenx, int leny,
                              int lox, int loy, int loz)
{
        int stride = blockDim.x*gridDim.x;

	for (int icell = blockDim.x*blockIdx.x+threadIdx.x;
	    icell < ncells; icell += stride) {
            int k =  icell /   (lenx*leny);
            int j = (icell - k*(lenx*leny)) /   lenx;
            int i = (icell - k*(lenx*leny)) - j*lenx;
            i += lox;
            j += loy;
            k += loz;
            //printf(" MyLaunchTest \n");
            //printf(" MyLaunchTest %d %d %d  \n", i,j,k);
            //printf(" MyLaunchTest %14.6e %14.6e \n", rho(i,j,k), temp(i,j,k));
	    gpu_RTY2W(i, j, k, rho, temp, mf, cdots);
        }
}


#endif
