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
      int oxy_idx = OXY_ID;
      extern_init(&(probin_file_name[0]),&probin_file_length,&fuel_idx,&oxy_idx);
    
      std::vector<int> npts(3,1);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 256;
      }
    
      Box domain(IntVect(D_DECL(0,0,0)),
                 IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

      std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
      for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
      }
    
      int max_size = 32;
      pp.query("max_size",max_size);
      BoxArray ba(domain);
      ba.maxSize(max_size);

      int num_spec;
      num_spec = NUM_SPECIES;

      DistributionMapping dm{ba};

      int num_grow = 0;
      MultiFab mass_frac(ba,dm,num_spec,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      //int cout_box = 0;
      for (MFIter mfi(mass_frac,tilesize); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			BL_TO_FORTRAN_N_3D(density[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
	//cout_box += 1;
      }
      //std::cout << "boxes ? " <<cout_box<<std::endl;

      ParmParse ppa("amr");
      std::string pltfile("plt");  
      ppa.query("plot_file",pltfile);
      std::string outfile = amrex::Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
      PlotFileFromMF(temperature,outfile);

      //Real time = 0.; pp.query("time",time);
      //Real dt=1.e-5; pp.query("dt",dt);

      MultiFab wdots(ba,dm,num_spec,num_grow);
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(mass_frac,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	const Box& box = mfi.tilebox();

	auto  mf      = mass_frac.array(mfi);
	auto  temp    = temperature.array(mfi);
	auto  rho     = density.array(mfi); 
	auto  cdots   = wdots.array(mfi);

	//FArrayBox& Y       = mass_frac[mfi];
	//FArrayBox& T       = temperature[mfi];
	//FArrayBox& D       = density[mfi];
	//FArrayBox& W       = wdots[mfi];

	//const auto len     = amrex::length(box);
	//const auto lo      = amrex::lbound(box);

	//const auto mf      = Y.view(lo); 
	//const auto temp    = T.view(lo); 
	//const auto rho     = D.view(lo); 
	//const auto cdots   = W.view(lo); 

	//for         (int k = 0; k < len.z; ++k) {
	//    for         (int j = 0; j < len.y; ++j) {
	//        for         (int i = 0; i < len.x; ++i) {
	//        
        //            eos.eos_RTY2W(rho(i,j,k), temp(i,j,k), mf(i,j,k,:), cdots(i,j,k,:));  
	//	}
	//    }
	//}

	amrex::ParallelFor(box,
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	    {
		gpu_RTY2W(i, j, k, rho, temp, mf, cdots);
		// put the gunction here
	    });

      }


      outfile = amrex::Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(wdots,outfile);

      extern_close();

    }

    Finalize();

    return 0;
}
