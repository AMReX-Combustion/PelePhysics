#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include "mechanism.h"
#include <GPU_misc.H>

#include <main_F.H>
#include <PlotFileFromMF.H>

#include "klu.h"


std::string inputs_name = "";

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {

      ParmParse pp;
    
      std::string probin_file = "probin";
      pp.query("probin_file",probin_file);
      int probin_file_length = probin_file.length();
      std::vector<int> probin_file_name(probin_file_length);

      for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

      int fuel_idx = FUEL_ID;
      int oxy_idx  = OXY_ID;
      int bath_idx = BATH_ID;

      extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx,&bath_idx);
    
      std::vector<int> npts(3,1);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 16;
      }
      npts[1] = 32;
    
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

      std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
      for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
      }
    
      int max_size = 16;
      pp.query("max_size",max_size);
      BoxArray ba(domain);
      ba.maxSize(max_size);

      int num_spec;
      num_spec = NUM_SPECIES;

      DistributionMapping dm{ba};

      int num_grow = 0;
      /* ORI QTY */
      MultiFab mass_frac(ba,dm,num_spec,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab energy(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);
      /* TMP QTY */
      MultiFab mass_frac_tmp(ba,dm,num_spec,num_grow);
      MultiFab temperature_tmp(ba,dm,1,num_grow);
      MultiFab energy_tmp(ba,dm,1,num_grow);
      MultiFab density_tmp(ba,dm,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));
    
      int box_count =0;
      for (MFIter mfi(mass_frac,tilesize); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			BL_TO_FORTRAN_N_3D(energy[mfi],0),
			BL_TO_FORTRAN_N_3D(density[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
	box_count +=1;
      }
      printf("That many boxes: %d \n", box_count);

      ParmParse ppa("amr");
      //std::string pltfile("plt");  
      //ppa.query("plot_file",pltfile);
      //std::string outfile = amrex::Concatenate(pltfile,0); // Need a number other than zero for reg test to pass

      std::string pltfile0("MF");  
      std::string outfile = amrex::Concatenate(pltfile0,0); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mass_frac,outfile);

      std::string pltfile1("TEMP");  
      outfile = amrex::Concatenate(pltfile1,0); // Need a number other than zero for reg test to pass
      PlotFileFromMF(temperature,outfile);

      Real time = 0.; pp.query("time",time);
      Real dt=1.e-07; pp.query("dt",dt);
      Real ndt=1000; pp.query("ndt",dt);
      MultiFab delta_t(ba,dm,1,num_grow);
      delta_t.setVal(dt,0,1,num_grow);

      MultiFab wdots(ba,dm,num_spec+1,num_grow);
      MultiFab systJ(ba,dm,(num_spec+1)*(num_spec+1),num_grow);
      MultiFab systRHS(ba,dm,num_spec+1,num_grow);
      MultiFab systRESNL(ba,dm,num_spec+1,num_grow);
      //MultiFab systNU(ba,dm,num_spec+1,num_grow);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	const Box& box = mfi.tilebox();
	int ncells = box.numPts();
        printf("Nb cells in box %d \n", ncells);

	/* dt */
	auto  deltas  = delta_t.array(mfi);
	/* Initial quantities (time n = q_n) */
	auto  mf        = mass_frac.array(mfi);
	auto  temp      = temperature.array(mfi);
	auto  nrgy      = energy.array(mfi);
	auto  rho       = density.array(mfi); 
	/* Tmp quantities (ite k = q_k) */
	auto  mf_tmp    = mass_frac_tmp.array(mfi);
	auto  temp_tmp  = temperature_tmp.array(mfi);
	auto  nrgy_tmp  = energy_tmp.array(mfi);
	auto  rho_tmp   = density_tmp.array(mfi); 
	/* RHS = q_k - cdots(q_k)*dt */
	auto  rhs       = systRHS.array(mfi);
	auto  cdots     = wdots.array(mfi);
	/* Jac of system = I-dt*J */
	auto  sJ        = systJ.array(mfi);
	/* nl residual for Newton iterations = RHS - q_n */
	auto  res_nl    = systRESNL.array(mfi);
	//auto  newton_update    = systNU.array(mfi);

        /* OBJECTS */
        klu_common *Common;
        klu_symbolic **Symbolic;
        klu_numeric **Numeric;
        Common   = new klu_common[ncells];
        Symbolic = new klu_symbolic*[ncells];
        Numeric  = new klu_numeric*[ncells];
        /* symbolic analysis */
        /* DENSE MAT right now */
        int* csr_col_count;
        int* csr_row_index;
	int* indx;
        amrex::Real* csr_val;
        amrex::Real* csr_b;
        amrex::Real* csr_x;
        int HP = 0;
	int offset_beg = 0;

        csr_col_count = (int *) malloc( (num_spec+2) * sizeof(int) );
        csr_row_index = (int *) malloc( (num_spec+1)*(num_spec+1)*sizeof(int) );
        indx          = (int *) malloc( (num_spec+1)*(num_spec+1)*sizeof(int) );
        csr_val       = (amrex::Real *) malloc( (num_spec+1)*(num_spec+1)*(ncells)*sizeof(amrex::Real) );
        csr_b         = (amrex::Real *) malloc( (num_spec+1)*(ncells)*sizeof(amrex::Real) );
        csr_x         = (amrex::Real *) malloc( (num_spec+1)*(ncells)*sizeof(amrex::Real) );

        SPARSITY_PREPROC_PRECOND(csr_row_index, csr_col_count, indx, &HP);

        for(int i = 0; i < ncells; ++i) { 
		// OPTI 1
		amrex::Real * csr_val_cell = csr_val + i * (num_spec+1) * (num_spec+1);
		// OPTI 2
		//amrex::Real csr_val_cell[(num_spec+1) * (num_spec+1)];
		//offset_beg = i * (num_spec+1) * (num_spec+1);
		//std::memcpy(csr_val_cell, csr_val + offset_beg, (num_spec+1) * (num_spec+1)*sizeof(amrex::Real));

		klu_defaults(&Common[i]);
		Symbolic[i]  = klu_analyze(num_spec+1, csr_col_count, csr_row_index, &Common[i]) ;
		Numeric[i]   = klu_factor(csr_col_count, csr_row_index, csr_val_cell, Symbolic[i], &Common[i]); 
  	}

        for (int stp = 0; stp < ndt; stp++) {
	        /* Copy init guess into q_k = q_0 */
	        amrex::ParallelFor(box,
	            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	            {
	        	gpu_CopyORI2TMP(i, j, k, rho, temp, nrgy, mf, rho_tmp, temp_tmp, nrgy_tmp, mf_tmp);
	            });

	        /* Estimate wdot for q_0 */
	        amrex::ParallelFor(box,
	            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	            {
	        	gpu_RTY2W(i, j, k, rho_tmp, temp_tmp, nrgy_tmp, mf_tmp, cdots);
	            });

	        /* Estimate RHS for q_0 */
	        amrex::ParallelFor(box,
	            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	            {
	        	gpu_RHS(i, j, k, rho_tmp, temp_tmp, mf_tmp, cdots, deltas, rhs);
	            });

	        /* Compute initial nl residual */
	        const auto lo  = amrex::lbound(box);
		const auto hi  = amrex::ubound(box);
	        const auto len = amrex::length(box);
		int icell = 0;
		for (int k = lo.z; k <= hi.z; ++k) {
			for (int j = lo.y; j <= hi.y; ++j) {
				AMREX_PRAGMA_SIMD
				for (int i = lo.x; i <= hi.x; ++i) {
					gpu_NLRES(i,j,k,icell, rho, temp, mf, rhs, res_nl, csr_b);
					icell = icell + 1;
				}
			}
		}

	        /* TODO COMPUTE NORM OF RES */

	        /* END INIT */

	        /* START LOOP */
	        bool newton_solved = false;
	        int newton_ite = 0;
	        //while (!newton_solved) {
	        while (newton_ite < 10) {
                        //printf("Ite number %d \n", newton_ite);
	        	newton_ite += 1;
	                /* Compute initial newton_update (delta q_k+1) */
			icell = 0;
			for (int k = lo.z; k <= hi.z; ++k) {
			for (int j = lo.y; j <= hi.y; ++j) {
			    AMREX_PRAGMA_SIMD
	                    for (int i = lo.x; i <= hi.x; ++i){
	                	gpu_resetNU(i, j, k, icell, csr_b, csr_x);
				icell = icell + 1;
	                    }
			}
			}

	                /* Jac chemistry */
	                amrex::ParallelFor(box,
	                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	                    {
	                	gpu_JAC(i, j, k, rho_tmp, temp_tmp, mf_tmp, sJ);
	                    });

	                /* Jac System */
			icell = 0;
			for (int k = lo.z; k <= hi.z; ++k) {
			for (int j = lo.y; j <= hi.y; ++j) {
			    AMREX_PRAGMA_SIMD
	                    for (int i = lo.x; i <= hi.x; ++i){
	                	gpu_J2SYSJ(i, j, k, icell, deltas, sJ, csr_val);
				icell = icell + 1;
	                    }
			}
			}

	                /* SOLVE LINEAR SYSTEM */
                        for(int i = 0; i < ncells; ++i) { 
	                	amrex::Real * csr_val_cell = csr_val + i * (num_spec+1) * (num_spec+1);
		                //amrex::Real csr_val_cell[(num_spec+1) * (num_spec+1)];
		                //offset_beg = i * (num_spec+1) * (num_spec+1);
		                //std::memcpy(csr_val_cell, csr_val + offset_beg, (num_spec+1) * (num_spec+1)*sizeof(amrex::Real));
                                amrex::Real * csr_x_cell = csr_x + i * (num_spec+1);
				//amrex::Real csr_x_cell[(num_spec+1)];
				//offset_beg = i * (num_spec+1);
				//std::memcpy(csr_x_cell, csr_x + offset_beg, (num_spec+1)*sizeof(amrex::Real));
	
	                	//klu_refactor(csr_col_count, csr_row_index, csr_val_cell, Symbolic[i], Numeric[i], &Common[i]);
		                Numeric[i]   = klu_factor(csr_col_count, csr_row_index, csr_val_cell, Symbolic[i], &Common[i]); 
				klu_solve(Symbolic[i], Numeric[i], num_spec+1, 1, csr_x_cell, &Common[i]) ; 
				//std::memcpy(csr_x + offset_beg, csr_x_cell, (num_spec+1)*sizeof(amrex::Real));
  	                }

	                /* UPDATE SOLUTION update q_tmp, it becomes q_(k+1, m+1)*/
			icell = 0;
			for (int k = lo.z; k <= hi.z; ++k) {
			for (int j = lo.y; j <= hi.y; ++j) {
			    AMREX_PRAGMA_SIMD
	                    for (int i = lo.x; i <= hi.x; ++i){
	                	gpu_UPDATETMP(i, j, k, icell, rho_tmp, temp_tmp, mf_tmp, csr_x);
				icell = icell + 1;
	                    }
			}
			}

	                /* Estimate wdot for q_k */
	                amrex::ParallelFor(box,
	                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	                    {
	        	        gpu_RTY2W(i, j, k, rho_tmp, temp_tmp, nrgy_tmp, mf_tmp, cdots);
	                    });

	                /* Estimate RHS for q_k */
	                amrex::ParallelFor(box,
	                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	                    {
	                	gpu_RHS(i, j, k, rho_tmp, temp_tmp, mf_tmp, cdots, deltas, rhs);
	                    });

	                /* Compute Newton nl residual */
			icell = 0;
			for (int k = lo.z; k <= hi.z; ++k) {
			for (int j = lo.y; j <= hi.y; ++j) {
			    AMREX_PRAGMA_SIMD
	                    for (int i = lo.x; i <= hi.x; ++i){
	                	gpu_NLRES(i, j, k, icell, rho, temp, mf, rhs, res_nl, csr_b);
				icell = icell + 1;
	                    }
			}
			}

	                ///* COMPUTE NORM OF RES */

	                ///* SEE IF WE CAN GET OUT OF NL SOLVE */
	                ///* BREAK LOOP */

	                ///* END LOOP */
	        } //( not newton_solved );

	        /* Copy q_tmp into q_(k+1) = iterations successful or other out criterion */
	        amrex::ParallelFor(box,
	            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	            {
	        	gpu_CopyTMP2ORI(i, j, k, rho_tmp, temp_tmp, nrgy_tmp, mf_tmp, rho, temp, nrgy, mf);
	            });
        }

      }

      //std::string pltfile0("MF");  
      outfile = amrex::Concatenate(pltfile0,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mass_frac,outfile);

      //std::string pltfile1("TEMP");  
      outfile = amrex::Concatenate(pltfile1,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(temperature,outfile);

      //std::string pltfile1("CDOTS");  
      //outfile = amrex::Concatenate(pltfile11,1); // Need a number other than zero for reg test to pass
      //PlotFileFromMF(systRHS,outfile);

      //std::string pltfile2("RESNL");  
      //outfile = amrex::Concatenate(pltfile2,2); // Need a number other than zero for reg test to pass
      //PlotFileFromMF(systRESNL,outfile);

      //std::string pltfile2("RESNL");  
      //outfile = amrex::Concatenate(pltfile2,2); // Need a number other than zero for reg test to pass
      //PlotFileFromMF(systRESNL,outfile);

      extern_close();

    }

    Finalize();

    return 0;
}
