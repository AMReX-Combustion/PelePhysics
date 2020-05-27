#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <main_F.H>
#include <PlotFileFromMF.H>
#include "mechanism.h"
#include <GPU_misc.H>

#include <EOS.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);
    {

      ParmParse pp;

      EOS::init();
    
      std::vector<int> npts(3,1);
      for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 128;
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

      ParmParse ppa("amr");
      std::string pltfile("plt");
      ppa.query("plot_file",pltfile);

      DistributionMapping dm{ba};

      int num_grow = 0;
      MultiFab mass_frac(ba,dm,NUM_SPECIES,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);
      MultiFab energy(ba,dm,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));
    
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mass_frac,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			BL_TO_FORTRAN_N_3D(density[mfi],0),
			BL_TO_FORTRAN_N_3D(energy[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));

	//const Box& gbox = mfi.tilebox();

	//Array4<Real> const& Y_a    = mass_frac.array(mfi);
	//Array4<Real> const& T_a    = temperature.array(mfi);
	//Array4<Real> const& rho_a  = density.array(mfi);
	//Array4<Real> const& e_a    = energy.array(mfi);

	//amrex::ParallelFor(gbox, 
	//    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	//	initialize_data(i, j, k, 
	//			Y_a, T_a, rho_a, e_a 
	//			dx, plo, phi);
	//});
      }

      MultiFab VarPlt(ba,dm,NUM_SPECIES+3,num_grow);
      MultiFab::Copy(VarPlt,mass_frac,0,0,NUM_SPECIES,num_grow);
      MultiFab::Copy(VarPlt,temperature,0,NUM_SPECIES,1,num_grow);
      MultiFab::Copy(VarPlt,density,0,NUM_SPECIES+1,1,num_grow);
      MultiFab::Copy(VarPlt,energy,0,NUM_SPECIES+2,1,num_grow);
      std::string initfile = amrex::Concatenate(pltfile,99); // Need a number other than zero for reg test to pass
      PlotFileFromMF(VarPlt,initfile);

      MultiFab cp(ba,dm,1,num_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mass_frac,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	const Box& gbox = mfi.tilebox();

	Array4<Real> const& Y_a    = mass_frac.array(mfi);
	Array4<Real> const& T_a    = temperature.array(mfi);
	Array4<Real> const& cp_a   = cp.array(mfi);

	amrex::ParallelFor(gbox, 
	    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
	    get_cp(i, j, k, 
	           Y_a, T_a, cp_a);
	});
      }

      EOS::close();

      std::string outfile = amrex::Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(cp,outfile);

    }

    amrex::Finalize();

    return 0;
}
