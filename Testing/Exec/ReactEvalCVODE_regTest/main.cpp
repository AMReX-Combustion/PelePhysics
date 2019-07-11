#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include <Transport_F.H>
#include <main_F.H>
#include <PlotFileFromMF.H>
#include <actual_Creactor.h>
#include "mechanism.h"

/**********************************/
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {

    int max_grid_size = 32;
    std::string probin_file="probin";
    std::string pltfile("plt");
    std::string txtfile_in=""; 
    /* CVODE inputs */
    int cvode_iE = 1;
    int ndt = 1; 
    Real dt = 1.e-5; 

    {
      /* ParmParse from the inputs file */
      ParmParse pp;
      
      // probin file
      pp.query("probin_file",probin_file);

      // domain size
      //pp.query("max_grid_size",max_grid_size);

      // final time
      pp.query("dt",dt);

      // time stepping
      pp.query("ndt",ndt); 

      pp.get("cvode_iE",cvode_iE);
      // Select CVODE type of energy employed.
      //1 for UV, 2 for HP
      //   1 = Internal energy
      //   anything else = enthalpy (PeleLM restart)
      
      pp.query("txtfile_in",txtfile_in);

    }

    amrex::Print() << "Integration method: ";
        amrex::Print() << "BDF (stiff)";
    amrex::Print() << std::endl;

    amrex::Print() << "Integration iteration method: ";
        amrex::Print() << "Newton";
    amrex::Print() << std::endl;
    
    amrex::Print() << std::endl;


    /* take care of probin init to initialize problem */
    int probin_file_length = probin_file.length();
    std::vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];

    extern_init(&(probin_file_name[0]),&probin_file_length, &cvode_iE);

#ifdef _OPENMP
#pragma omp parallel
#endif  
    /* Initialize D/CVODE reactor */
    int cvode_ncells = 1;
    reactor_init(&cvode_iE, &cvode_ncells);

    /* make domain and BoxArray */
    std::vector<int> npts(3,1);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 128;
    }

    Box domain(IntVect(D_DECL(0,0,0)),
	       IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    /* Additional defs to initialize domain */
    std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
    for (int i=0; i<BL_SPACEDIM; ++i) {
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
    }

    int Ncomp;
    Ncomp =  NUM_SPECIES;

    /* Create MultiFabs with no ghost cells */
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, Ncomp+1, 0);
    MultiFab rY_source_ext(ba,dm,Ncomp,0);
    MultiFab mfE(ba, dm, 1, 0);
    MultiFab rY_source_energy_ext(ba,dm,1,0);
    MultiFab temperature(ba,dm,1,0);

    IntVect tilesize(D_DECL(10240,8,32));

    /* INITIALIZE DATA */
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi ){
        const Box& box = mfi.tilebox();
        initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
               		BL_TO_FORTRAN_N_3D(mf[mfi],0),
               		BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
		        BL_TO_FORTRAN_N_3D(mfE[mfi],0),
		        BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
    }
     
    ParmParse ppa("amr");
    ppa.query("plot_file",pltfile);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf,tilesize); mfi.isValid(); ++mfi )
    {
        const Box& box = mfi.tilebox();

	FArrayBox& Fb     = mf[mfi];
	FArrayBox& Fbsrc  = rY_source_ext[mfi];
	FArrayBox& FbE    = mfE[mfi];
	FArrayBox& FbEsrc = rY_source_energy_ext[mfi];

        /* ADVANCE */
        int reInit = 1;
        // NOT USED
        double pressure = 1013250.0;

        /* Pack the data */
	// rhoY,T
	double tmp_vect[Ncomp+1];
	// rhoY_src_ext
	double tmp_src_vect[Ncomp];
	// rhoE/rhoH
	double tmp_vect_energy[1];
	double tmp_src_vect_energy[1];

	for (BoxIterator bit(box); bit.ok(); ++bit) {
		/* Fill the vectors */
		tmp_vect_energy[0]     = FbE(bit(),0);
		tmp_src_vect_energy[0] = FbEsrc(bit(),0);
		for (int i=0;i<Ncomp; i++){
			tmp_vect[i]     = Fb(bit(),i);
			tmp_src_vect[i] = Fbsrc(bit(),i);
		}
	        tmp_vect[Ncomp] = Fb(bit(), Ncomp);
                /* Solve the problem */
	        Real time_tmp, dt_incr;
	        dt_incr =  dt / ndt;
	        time_tmp = 0.0;
	        for (int i = 0; i < ndt; ++i) {
		        react(tmp_vect, tmp_src_vect, 
				tmp_vect_energy, tmp_src_vect_energy,
				&pressure, &dt_incr, &time_tmp,
				&reInit);
		        // fix new dt_incr to chosen value, hoping cvode will reach it
		        dt_incr = dt / ndt;
	        }

                /* Unpack the data ? */
		for (int i=0;i<Ncomp+1; i++){
			Fb(bit(),i) = tmp_vect[i];
		}
		FbE(bit(),0) = tmp_vect_energy[0];
	}
    }



    std::string outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
    MultiFab::Copy(temperature,mf,Ncomp,0,1,0);
    //PlotFileFromMF(mf,outfile);
    PlotFileFromMF(temperature,outfile);

    extern_close();

    }

    amrex::Finalize();

    return 0;
}
