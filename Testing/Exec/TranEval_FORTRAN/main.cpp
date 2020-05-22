#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <Transport_F.H>
#include <main_F.H>
#include <PlotFileFromMF.H>
#include "mechanism.h"

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

      extern_init(&(probin_file_name[0]),&probin_file_length);
    
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

      DistributionMapping dm{ba};

      int num_grow = 0;
      MultiFab mass_frac(ba,dm,NUM_SPECIES,num_grow);
      MultiFab temperature(ba,dm,1,num_grow);
      MultiFab density(ba,dm,1,num_grow);

      IntVect tilesize(D_DECL(10240,8,32));
    
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mass_frac,tilesize); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			BL_TO_FORTRAN_N_3D(density[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
      }

      MultiFab D(ba,dm,NUM_SPECIES+3,num_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mass_frac,tilesize); mfi.isValid(); ++mfi) {
	const Box& box = mfi.tilebox();
	get_transport_coeffs_F(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
			     BL_TO_FORTRAN_N_3D(mass_frac[mfi],0),
			     BL_TO_FORTRAN_N_3D(temperature[mfi],0),
			     BL_TO_FORTRAN_N_3D(density[mfi],0),
			     BL_TO_FORTRAN_N_3D(D[mfi],0),
			     BL_TO_FORTRAN_N_3D(D[mfi],NUM_SPECIES),
			     BL_TO_FORTRAN_N_3D(D[mfi],NUM_SPECIES+1),
			     BL_TO_FORTRAN_N_3D(D[mfi],NUM_SPECIES+2));
      }

      ParmParse ppa("amr");
      std::string pltfile("plt");  ppa.query("plot_file",pltfile);
      std::string outfile = amrex::Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(D,outfile);

      extern_close();

    }

    amrex::Finalize();

    return 0;
}
