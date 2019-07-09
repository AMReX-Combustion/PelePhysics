#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <Fuego_EOS.H>
#include "chemistry_file.H"

AMREX_GPU_DEVICE
inline
void
gpu_RTY2W(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& nrgy,
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& cdots) noexcept
{
	amrex::Real mw[NUM_SPECIES];
	amrex::GpuArray<amrex::Real,NUM_SPECIES> mf_pt;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> ei_pt;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> cdots_pt;
	amrex::Real Cv_pt;

	/* Specs */
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
		cdots_pt[n] = cdots(i,j,k,n);
	}

	EOS eos; 
        eos.eos_EY2T(mf_pt.arr, nrgy(i,j,k), temp(i,j,k));
	eos.eos_RTY2W(rho(i,j,k), temp(i,j,k), mf_pt.arr, cdots_pt.arr);
	eos.eos_T2EI(temp(i,j,k), ei_pt.arr);
	eos.eos_TY2Cv(temp(i,j,k), mf_pt.arr, &Cv_pt);
	get_mw(mw);
	
	cdots(i,j,k,NUM_SPECIES) = 0.0;
	for (int n = 0; n < NUM_SPECIES; n++) {
		cdots(i,j,k,n) = cdots_pt[n] * mw[n];
		cdots(i,j,k,NUM_SPECIES) = cdots(i,j,k,NUM_SPECIES) - ei_pt[n] * cdots(i,j,k,n);
	}
	cdots(i,j,k,NUM_SPECIES) = cdots(i,j,k,NUM_SPECIES) / rho(i,j,k) / Cv_pt; 
}

AMREX_GPU_DEVICE
inline
void
gpu_JAC(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& sJ) noexcept
{
	amrex::Real mw[NUM_SPECIES];
	amrex::GpuArray<amrex::Real,NUM_SPECIES> mf_pt;
	amrex::GpuArray<amrex::Real,(NUM_SPECIES+1)*(NUM_SPECIES+1)> J_pt;

	/* Specs */
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_pt[n] = mf(i,j,k,n);
	}
	for (int n = 0; n < NUM_SPECIES+1; n++) {
		for (int m = 0; m < NUM_SPECIES+1; m++) {
			J_pt[(NUM_SPECIES+1)*n+m] = 0.0; //sJ(i,j,k,(NUM_SPECIES+1)*n+m);
		}
	}

	EOS eos; 
	eos.eos_RTY2JAC(rho(i,j,k), temp(i,j,k), mf_pt.arr, J_pt.arr);
	get_mw(mw);

	for (int n = 0; n < NUM_SPECIES; n++) {
		for (int m = 0; m < NUM_SPECIES; m++) {
			sJ(i,j,k,(NUM_SPECIES+1)*n+m) = J_pt[(NUM_SPECIES+1)*n+m]*mw[m]/mw[n];
		}
		sJ(i,j,k,(NUM_SPECIES+1)*n+NUM_SPECIES) = J_pt[(NUM_SPECIES+1)*n+NUM_SPECIES]/mw[n];
		sJ(i,j,k,(NUM_SPECIES+1)*NUM_SPECIES+n) = J_pt[(NUM_SPECIES+1)*NUM_SPECIES+n]*mw[n];
	}
	sJ(i,j,k,(NUM_SPECIES+1)*(NUM_SPECIES+1)-1) = J_pt[(NUM_SPECIES+1)*(NUM_SPECIES+1)-1];
}


AMREX_GPU_DEVICE
inline
void
gpu_J2SYSJ(int i, int j, int k,
                int icell,  
		amrex::Array4<amrex::Real> const& deltat,
		amrex::Array4<amrex::Real> const& sJ,
                amrex::Real * csr_val) noexcept
{
        amrex::Real * csr_val_cell = csr_val + icell * (NUM_SPECIES+1)*(NUM_SPECIES+1);
	for (int n = 0; n < NUM_SPECIES+1; n++) {
		sJ(i,j,k,(NUM_SPECIES+1)*n+n) = 1 - deltat(i,j,k) * sJ(i,j,k,(NUM_SPECIES+1)*n+n);
                csr_val_cell[(NUM_SPECIES+1)*n+n] = sJ(i,j,k,(NUM_SPECIES+1)*n+n);
		for (int m = 0; m < n; m++) { 
			sJ(i,j,k,(NUM_SPECIES+1)*n+m) = - deltat(i,j,k) * sJ(i,j,k,(NUM_SPECIES+1)*n+m);
                        csr_val_cell[(NUM_SPECIES+1)*m+n] = sJ(i,j,k,(NUM_SPECIES+1)*n+m);
		}	
		for (int m = n+1; m < (NUM_SPECIES+1); m++) { 
			sJ(i,j,k,(NUM_SPECIES+1)*n+m) = - deltat(i,j,k) * sJ(i,j,k,(NUM_SPECIES+1)*n+m);
                        csr_val_cell[(NUM_SPECIES+1)*m+n] = sJ(i,j,k,(NUM_SPECIES+1)*n+m);
		}	
	}
        //printf(" J 3rd line ? \n");
        //for (int n = 0; n < NUM_SPECIES+1; n++) {
        //    printf(" %d %14.6e  %14.6e\n", n,  csr_val_cell[2*(NUM_SPECIES+1) + n], sJ(i,j,k,(NUM_SPECIES+1)*n + 2));
        //}
}

AMREX_GPU_DEVICE
inline
void
gpu_RHS(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& cdots,
		amrex::Array4<amrex::Real> const& deltat,
		amrex::Array4<amrex::Real> const& rhs) noexcept
{

	/* RHS is y_k */
	for (int n = 0; n < NUM_SPECIES; n++) {
		rhs(i,j,k,n) = rho(i,j,k) * mf(i,j,k,n) - deltat(i,j,k) * cdots(i,j,k,n);
	}
	rhs(i,j,k,NUM_SPECIES) = temp(i,j,k) - deltat(i,j,k) * cdots(i,j,k,NUM_SPECIES);
}

AMREX_GPU_DEVICE
inline
void
gpu_NLRES(int i, int j, int k,  
                int icell,
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& rhs,
		amrex::Array4<amrex::Real> const& res_nl,
                amrex::Real * csr_b) noexcept
{
        amrex::Real * csr_b_cell = csr_b + icell * (NUM_SPECIES+1);
	/* RHS is y_k */
	for (int n = 0; n < NUM_SPECIES; n++) {
		res_nl(i,j,k,n) = -(rhs(i,j,k,n) - rho(i,j,k) * mf(i,j,k,n));
                csr_b_cell[n] = res_nl(i,j,k,n);
	}
	res_nl(i,j,k,NUM_SPECIES) = -(rhs(i,j,k,NUM_SPECIES) - temp(i,j,k));
        csr_b_cell[NUM_SPECIES] = res_nl(i,j,k,NUM_SPECIES); 
}

AMREX_GPU_DEVICE
inline
void
gpu_UPDATETMP(int i, int j, int k,  
                int icell,
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& mf,
                amrex::Real * csr_x) noexcept
{
        amrex::Real * csr_x_cell = csr_x + icell * (NUM_SPECIES+1);
	/* RHS is y_k */
        amrex::Real rho_tmp = 0.0;
	for (int n = 0; n < NUM_SPECIES; n++) {
                /* mf temporarily contains rho Y */
		mf(i,j,k,n) = rho(i,j,k)*mf(i,j,k,n) + csr_x_cell[n];
                rho_tmp  = rho_tmp + mf(i,j,k,n); 
	}
        rho(i,j,k) = rho_tmp;
	for (int n = 0; n < NUM_SPECIES; n++) {
                /* mf temporarily contains rho Y */
		mf(i,j,k,n) = mf(i,j,k,n)/rho(i,j,k);
	}
	temp(i,j,k) = temp(i,j,k) + csr_x_cell[NUM_SPECIES];
        //printf("Update T and T and rho %14.6e %14.6e %14.6e \n", csr_x_cell[NUM_SPECIES], temp(i,j,k), rho(i,j,k));
}


AMREX_GPU_DEVICE
inline
void
gpu_CopyTMP2ORI(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rho_t,
		amrex::Array4<amrex::Real> const& temp_t,
		amrex::Array4<amrex::Real> const& nrgy_t,
		amrex::Array4<amrex::Real> const& mf_t,
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& nrgy,
		amrex::Array4<amrex::Real> const& mf) noexcept
{
	/* RHS is y_k */
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf(i,j,k,n) = mf_t(i,j,k,n);
	}
	temp(i,j,k) = temp_t(i,j,k);
	nrgy(i,j,k) = nrgy_t(i,j,k);
	rho(i,j,k)  = rho_t(i,j,k);
        //printf("i,j,k, T, rho, E %d %d %d %14.6e %14.6e %14.6e\n", i, j, k, temp(i,j,k), rho(i,j,k), nrgy(i,j,k));
}

AMREX_GPU_DEVICE
inline
void
gpu_CopyORI2TMP(int i, int j, int k,  
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& nrgy,
		amrex::Array4<amrex::Real> const& mf,
		amrex::Array4<amrex::Real> const& rho_t,
		amrex::Array4<amrex::Real> const& temp_t,
		amrex::Array4<amrex::Real> const& nrgy_t,
		amrex::Array4<amrex::Real> const& mf_t) noexcept
{
	/* RHS is y_k */
	for (int n = 0; n < NUM_SPECIES; n++) {
		mf_t(i,j,k,n) = mf(i,j,k,n);
	}
	temp_t(i,j,k) = temp(i,j,k);
	nrgy_t(i,j,k) = nrgy(i,j,k);
	rho_t(i,j,k)  = rho(i,j,k);
}

AMREX_GPU_DEVICE
inline
void
gpu_resetNU(int i, int j, int k,  
                int icell,
                amrex::Real * csr_x) noexcept
		//amrex::Array4<amrex::Real> const& nu) noexcept
{
        amrex::Real * csr_x_cell = csr_x + icell * (NUM_SPECIES+1);
	for (int n = 0; n < NUM_SPECIES+1; n++) {
	//	nu(i,j,k,n) = 0.0;
                csr_x_cell[n] = 0.0;
	}
}

#endif
