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
      std::string pltfile("plt");  
      ppa.query("plot_file",pltfile);
      std::string outfile = amrex::Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mass_frac,outfile);

      std::string pltfile1("TEMP");  
      outfile = amrex::Concatenate(pltfile1,0); // Need a number other than zero for reg test to pass
      PlotFileFromMF(temperature,outfile);

      Real time = 0.; pp.query("time",time);
      Real dt=1.e-07; pp.query("dt",dt);
      Real ndt=500; pp.query("ndt",dt);
      //Real ndt=1000; pp.query("ndt",dt);
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
        printf("Nb cells in box %d \n", box.numPts());

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

        /* cuSolver */
        /* FLAGS */
        /* OBJECTS */
        klu_common *Common;
        klu_symbolic **Symbolic;
        klu_numeric **Numeric;
        Common = = new klu_common[box.numPts()];
        Symbolic = new klu_symbolic*[box.numPts()];
        Numeric = new klu_numeric*[box.numPts()];

        for(int i = 0; i < NCELLS; ++i) { 
}

        /* symbolic analysis */
        /* DENSE MAT right now */
        int* csr_col_count;
        int* csr_row_index;
        amrex::Real* csr_val;
        amrex::Real* csr_b;
        amrex::Real* csr_x;
        int HP = 0;

        csr_col_count = (int *) malloc( (num_spec+2) * sizeof(int) );
        csr_row_index = (int *) malloc( (num_spec+1)*(num_spec+1)*sizeof(int) );
        csr_val       = (amrex::Real *) malloc( (num_spec+1)*(num_spec+1)*(box.numPts())*sizeof(amrex::Real) );
        csr_b         = (amrex::Real *) malloc( (num_spec+1)*(box.numPts())*sizeof(amrex::Real) );
        csr_x         = (amrex::Real *) malloc( (num_spec+1)*(box.numPts())*sizeof(amrex::Real) );

        SPARSITY_PREPROC_PRECOND(csr_col_count, csr_row_index, &HP);
        //for  (int i = 0; i < num_spec+2; i++) {
        //  printf(" csr_row_count %d %d \n", i, csr_row_count[i]);
        //}

        cusolver_status = cusolverSpXcsrqrAnalysisBatched(cusolverHandle,
                                                        num_spec+1,
                                                        num_spec+1,
                                                        (num_spec+1)*(num_spec+1),
                                                        descrA,
                                                        csr_row_count,
                                                        csr_col_index,
                                                        info);
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

        /* allocate working space */
        cusolver_status = cusolverSpDcsrqrBufferInfoBatched(cusolverHandle,
                                                        num_spec+1,
                                                        num_spec+1,
                                                        (num_spec+1)*(num_spec+1),
                                                        descrA,
                                                        csr_val,
                                                        csr_row_count,
                                                        csr_col_index,
                                                        box.numPts(),
                                                        info,
                                                        &internalDataInBytes,
                                                        &workspaceInBytes);
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

        cudaStat1 = cudaMalloc((void**)&buffer_qr, workspaceInBytes);
        assert(cudaStat1 == cudaSuccess);

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
	        int ncells = box.numPts();
	        const auto lo  = amrex::lbound(box);
	        const auto len = amrex::length(box);
	        const auto ec = Gpu::ExecutionConfig(ncells);
	        amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	        [=] AMREX_GPU_DEVICE () noexcept {
	            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	                icell < ncells; icell += stride) {
	                int k =  icell /   (len.x*len.y);
	        	int j = (icell - k*(len.x*len.y)) /   len.x;
	        	int i = (icell - k*(len.x*len.y)) - j*len.x;
	        	i += lo.x;
	        	j += lo.y;
	        	k += lo.z;
	        	gpu_NLRES(i, j, k, icell, rho, temp, mf, rhs, res_nl, csr_b);
	            }
	        });

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
	                amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	                [=] AMREX_GPU_DEVICE () noexcept {
	                    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	                        icell < ncells; icell += stride) {
	                        int k =  icell /   (len.x*len.y);
	                	int j = (icell - k*(len.x*len.y)) /   len.x;
	                	int i = (icell - k*(len.x*len.y)) - j*len.x;
	                	i += lo.x;
	                	j += lo.y;
	                	k += lo.z;
	                	gpu_resetNU(i, j, k, icell, csr_x);
	                    }
	                });

	                /* Jac chemistry */
	                amrex::ParallelFor(box,
	                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	                    {
	                	gpu_JAC(i, j, k, rho_tmp, temp_tmp, mf_tmp, sJ);
	                    });

	                /* Jac System */
	                amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	                [=] AMREX_GPU_DEVICE () noexcept {
	                    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	                        icell < ncells; icell += stride) {
	                        int k =  icell /   (len.x*len.y);
	                	int j = (icell - k*(len.x*len.y)) /   len.x;
	                	int i = (icell - k*(len.x*len.y)) - j*len.x;
	                	i += lo.x;
	                	j += lo.y;
	                	k += lo.z;
	                	gpu_J2SYSJ(i, j, k, icell, deltas, sJ, csr_val);
	                    }
	                });

	                /* SOLVE LINEAR SYSTEM */
                        /* allocate working space */
                        cusolver_status = cusolverSpDcsrqrBufferInfoBatched(cusolverHandle,
                                                                        num_spec+1,
                                                                        num_spec+1,
                                                                        (num_spec+1)*(num_spec+1),
                                                                        descrA,
                                                                        csr_val,
                                                                        csr_row_count,
                                                                        csr_col_index,
                                                                        box.numPts(),
                                                                        info,
                                                                        &internalDataInBytes,
                                                                        &workspaceInBytes);
                        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

                        //cudaStat1 = cudaMalloc((void**)&buffer_qr, workspaceInBytes);
                        //assert(cudaStat1 == cudaSuccess);

                        cusolver_status = cusolverSpDcsrqrsvBatched(cusolverHandle,
                                                                        num_spec+1,
                                                                        num_spec+1,
                                                                        (num_spec+1)*(num_spec+1),
                                                                        descrA,
                                                                        csr_val,
                                                                        csr_row_count,
                                                                        csr_col_index,
                                                                        csr_b,
                                                                        csr_x,
                                                                        box.numPts(),
                                                                        info,
                                                                        buffer_qr);

	                /* UPDATE SOLUTION update q_tmp, it becomes q_(k+1, m+1)*/
	                amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	                [=] AMREX_GPU_DEVICE () noexcept {
	                    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	                        icell < ncells; icell += stride) {
	                        int k =  icell /   (len.x*len.y);
	                	int j = (icell - k*(len.x*len.y)) /   len.x;
	                	int i = (icell - k*(len.x*len.y)) - j*len.x;
	                	i += lo.x;
	                	j += lo.y;
	                	k += lo.z;
	                	gpu_UPDATETMP(i, j, k, icell, rho_tmp, temp_tmp, mf_tmp, csr_x);
	                    }
	                });


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
	                amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, amrex::Gpu::gpuStream()>>>(
	                [=] AMREX_GPU_DEVICE () noexcept {
	                    for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	                        icell < ncells; icell += stride) {
	                        int k =  icell /   (len.x*len.y);
	                	int j = (icell - k*(len.x*len.y)) /   len.x;
	                	int i = (icell - k*(len.x*len.y)) - j*len.x;
	                	i += lo.x;
	                	j += lo.y;
	                	k += lo.z;
	                	gpu_NLRES(i, j, k, icell, rho, temp, mf, rhs, res_nl, csr_b);
	                    }
	                });

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

      std::string pltfile0("MF");  
      outfile = amrex::Concatenate(pltfile0,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mass_frac,outfile);

      //std::string pltfile1("TEMP");  
      outfile = amrex::Concatenate(pltfile1,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(temperature,outfile);

      //std::string pltfile11("CDOTS");  
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
