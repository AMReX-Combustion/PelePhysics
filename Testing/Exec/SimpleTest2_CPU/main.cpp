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

///* CUDA cuSolver */
//#include <cuda_runtime.h>
//#include <cublas_v2.h>
//#include <cusolverSp.h>
//#include <cusparse.h> 
//#include <cusparse.h> 

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
	npts[i] = 8;
      }
      //npts[1] = 16;
    
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

      std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
      for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
      }
    
      int max_size = 4;
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
      MultiFab density(ba,dm,1,num_grow);
      /* TMP QTY */
      MultiFab mass_frac_tmp(ba,dm,num_spec,num_grow);
      MultiFab temperature_tmp(ba,dm,1,num_grow);
      MultiFab density_tmp(ba,dm,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));
    
      int box_count =0;
      for (MFIter mfi(mass_frac,tilesize); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
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

      Real time = 0.; pp.query("time",time);
      Real dt=1.e-10; pp.query("dt",dt);
      //Real dt=1.; pp.query("dt",dt);
      MultiFab delta_t(ba,dm,1,num_grow);
      delta_t.setVal(dt,0,1,num_grow);

      MultiFab wdots(ba,dm,num_spec+1,num_grow);
      MultiFab systJ(ba,dm,(num_spec+1)*(num_spec+1),num_grow);
      MultiFab systRHS(ba,dm,num_spec+1,num_grow);
      MultiFab systRESNL(ba,dm,num_spec+1,num_grow);
      MultiFab systNU(ba,dm,num_spec+1,num_grow);

    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	const Box& box = mfi.tilebox();

	/* dt */
	auto  deltas  = delta_t.array(mfi);
	/* Initial quantities (time n = q_n) */
	auto  mf        = mass_frac.array(mfi);
	auto  temp      = temperature.array(mfi);
	auto  rho       = density.array(mfi); 
	/* Tmp quantities (ite k = q_k) */
	auto  mf_tmp    = mass_frac_tmp.array(mfi);
	auto  temp_tmp  = temperature_tmp.array(mfi);
	auto  rho_tmp   = density_tmp.array(mfi); 
	/* RHS = q_k - cdots(q_k)*dt */
	auto  rhs       = systRHS.array(mfi);
	auto  cdots     = wdots.array(mfi);
	/* Jac of system = I-dt*J */
	auto  sJ        = systJ.array(mfi);
	/* nl residual for Newton iterations = RHS - q_n */
	auto  res_nl    = systRESNL.array(mfi);
	auto  newton_update    = systNU.array(mfi);

	/* Copy init guess into q_k = q_0 */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_CopyORI2TMP(i, j, k, rho, temp, mf, rho_tmp, temp_tmp, mf_tmp);
	    });

	/* Estimate wdot for q_0 */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_RTY2W(i, j, k, rho_tmp, temp_tmp, mf_tmp, cdots);
	    });

	/* Estimate RHS for q_0 */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_RHS(i, j, k, rho_tmp, temp_tmp, mf_tmp, cdots, deltas, rhs);
	    });

	/* Compute initial nl residual */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_NLRES(i, j, k, rho, temp, mf, rhs, res_nl);
	    });

	/* COMPUTE NORM OF RES */
	/* END INIT */

	/* START LOOP */
	/* Compute initial newton_update (delta q_k+1) */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_resetNU(i, j, k, newton_update);
	    });

	/* Jac chemistry */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_JAC(i, j, k, rho_tmp, temp_tmp, mf_tmp, sJ);
	    });

	/* Jac System */
	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_J2SYSJ(i, j, k, deltas, sJ);
	    });

	/* SOLVE LINEAR SYSTEM */

	///* UPDATE SOLUTION update q_tmp, it becomes q_(k+1, m+1)*/


	///* Estimate wdot for q_k */
	//amrex::ParallelFor(box,
	//    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	//    {
	//	gpu_RTY2W(i, j, k, rho_tmp, temp_tmp, mf_tmp, cdots);
	//    });

	///* Estimate RHS for q_k */
	//amrex::ParallelFor(box,
	//    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	//    {
	//	gpu_RHS(i, j, k, rho_tmp, temp_tmp, mf_tmp, cdots, deltas, rhs);
	//    });

	///* Compute Newton nl residual */
	//amrex::ParallelFor(box,
	//    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	//    {
	//	gpu_NLRES(i, j, k, rho, temp, mf, rhs, res_nl);
	//    });

	///* COMPUTE NORM OF RES */

	///* SEE IF WE CAN GET OUT OF NL SOLVE */
	///* BREAK LOOP */

	///* END LOOP */

	///* Copy q_tmp into q_(k+1) = iterations successful or other out criterion */
	//amrex::ParallelFor(box,
	//    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	//    {
	//	gpu_CopyTMP2ORI(i, j, k, rho_tmp, temp_tmp, mf_tmp, rho, temp, mf);
	//    });
	//    });

	/* COMPUTE NORM OF RES */
	/* END INIT */

      }

      std::string pltfile0("MF");  
      outfile = amrex::Concatenate(pltfile0,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mass_frac,outfile);

      std::string pltfile1("CDOTS");  
      outfile = amrex::Concatenate(pltfile1,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(wdots,outfile);

      std::string pltfile11("RHS");  
      outfile = amrex::Concatenate(pltfile11,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(systRHS,outfile);

      std::string pltfile2("RESNL");  
      outfile = amrex::Concatenate(pltfile2,2); // Need a number other than zero for reg test to pass
      PlotFileFromMF(systRESNL,outfile);

      extern_close();

    }

    Finalize();

    return 0;
}
