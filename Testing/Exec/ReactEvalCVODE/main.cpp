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
#ifdef AMREX_USE_SUNDIALS_3x4x
#include <actual_Creactor.h>
#else
#include <actual_reactor.H> 
#endif

/**********************************/
int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {

    int max_grid_size = 16;
    std::string probin_file="probin";
    std::string pltfile("plt");
    /* CVODE inputs */
    int cvode_ncells = 1;
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
         
      // nb of cells to integrate
      pp.query("cvode_ncells",cvode_ncells);

    }

    amrex::Print() << "Integration method: ";
        amrex::Print() << "BDF (stiff)";
    amrex::Print() << std::endl;

    amrex::Print() << "Integration iteration method: ";
        amrex::Print() << "Newton";
    amrex::Print() << std::endl;

    amrex::Print() << "Type of reactor: ";
        amrex::Print() << cvode_iE;
    amrex::Print() << std::endl;

    amrex::Print() << "Integrating 2x1024x2 box for: ";
        amrex::Print() << dt << " seconds";
    amrex::Print() << std::endl;


    amrex::Print() << std::endl;

    /* take care of probin init to initialize problem */
    int probin_file_length = probin_file.length();
    std::vector<int> probin_file_name(probin_file_length);
    for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];
    extern_init(&(probin_file_name[0]),&probin_file_length, &cvode_iE);

    /* Initialize D/CVODE reactor */
    reactor_init(&cvode_iE, &cvode_ncells);

    /* make domain and BoxArray */
    std::vector<int> npts(3,1);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	npts[i] = 2;
    }
    npts[1] = 1024;

    Box domain(IntVect(D_DECL(0,0,0)),
	       IntVect(D_DECL(npts[0]-1,npts[1]-1,npts[2]-1)));

    BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    /* Additional defs to initialize domain */
    std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
    for (int i=0; i<BL_SPACEDIM; ++i) {
	plo[i] = 0.0; //(i+1)*0.35;
	phi[i] = domain.length(i);
	dx[i] = (phi[i] - plo[i])/domain.length(i);
    }

    int Ncomp;
    get_num_spec(&Ncomp);

    /* Create MultiFabs with no ghost cells */
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, Ncomp+1, 0);
    MultiFab rY_source_ext(ba,dm,Ncomp,0);
    MultiFab mfE(ba, dm, 1, 0);
    MultiFab rY_source_energy_ext(ba,dm,1,0);
    MultiFab temperature(ba,dm,1,0);
    MultiFab fctCount(ba,dm,1,0);

    /* INITIALIZE DATA */
#ifdef _OPENMP
#pragma omp parallel
#endif
    int count_mf = 0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi ){
        count_mf = count_mf + 1;	
        const Box& box = mfi.tilebox();
        initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
               		BL_TO_FORTRAN_N_3D(mf[mfi],0),
               		BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
		        BL_TO_FORTRAN_N_3D(mfE[mfi],0),
		        BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
			&(dx[0]), &(plo[0]), &(phi[0]));
    }
    amrex::Print() << "Number of boxes ? (64): \n" << count_mf;
    std::string outfile = Concatenate(pltfile,0); // Need a number other than zero for reg test to pass
    //// Specs
    PlotFileFromMF(mf,outfile);

     
    ParmParse ppa("amr");
    ppa.query("plot_file",pltfile);

    /* ADVANCE */
    Real time = 0.0;
    int reInit = 1;
    Real dt_incr   = dt/ndt;
    // not used anyway
    double pressure = 1013250.0;

    amrex::Print() << " \n STARTING THE ADVANCE \n";

    count_mf = 0;
    for ( MFIter mfi(mf); mfi.isValid(); ++mfi )
    {
	/* Prints to follow the computation */
        count_mf = count_mf + 1;	
        amrex::Print() << "Treating box " << count_mf << std::endl;

        const Box& box = mfi.tilebox();

	FArrayBox& Fb     = mf[mfi];
	FArrayBox& Fbsrc  = rY_source_ext[mfi];
	FArrayBox& FbE    = mfE[mfi];
	FArrayBox& FbEsrc = rY_source_energy_ext[mfi];
	FArrayBox& Fct    = fctCount[mfi];

	const auto len     = amrex::length(box);
	const auto lo      = amrex::lbound(box);

	const auto rhoY    = Fb.view(lo);
	const auto rhoE    = FbE.view(lo);
	const auto frcExt  = Fbsrc.view(lo); 
	const auto frcEExt = FbEsrc.view(lo);
	const auto fc      = Fct.view(lo); 

        /* Pack the data */
	// rhoY,T
	double tmp_vect[cvode_ncells*(Ncomp+1)];
	// rhoY_src_ext
	double tmp_src_vect[cvode_ncells*(Ncomp)];
	// rhoE/rhoH
	double tmp_vect_energy[cvode_ncells];
	double tmp_src_vect_energy[cvode_ncells];

	int indx_i[cvode_ncells];
	int indx_j[cvode_ncells];
	int indx_k[cvode_ncells];
	
	int nc = 0;
	int num_cell_cvode_int = 0;
	for         (int k = 0; k < len.z; ++k) {
	    for         (int j = 0; j < len.y; ++j) {
	        for         (int i = 0; i < len.x; ++i) {
		    /* Fill the vectors */
	            for (int sp=0;sp<Ncomp; sp++){
	                tmp_vect[nc*(Ncomp+1) + sp]   = rhoY(i,j,k,sp);
		        tmp_src_vect[nc*Ncomp + sp]   = frcExt(i,j,k,sp);
		    }
		    tmp_vect[nc*(Ncomp+1) + Ncomp]    = rhoY(i,j,k,Ncomp);
		    tmp_vect_energy[nc]               = rhoE(i,j,k,0);
		    tmp_src_vect_energy[nc]           = frcEExt(i,j,k,0);
		    //
		    indx_i[nc] = i;
		    indx_j[nc] = j;
		    indx_k[nc] = k;
		    //
		    nc = nc+1;
		    //
		    num_cell_cvode_int = num_cell_cvode_int + 1;
		    if (nc == cvode_ncells) {
			time = 0.0;
			dt_incr =  dt/ndt;
                        reInit = 1;
			for (int ii = 0; ii < ndt; ++ii) {
	                    fc(i,j,k) = react(tmp_vect, tmp_src_vect,
		                tmp_vect_energy, tmp_src_vect_energy,
		                &pressure, &dt_incr, &time,
				&reInit);
		            dt_incr =  dt/ndt;
			    reInit = 1;
			}
			//printf(" time reached %14.6e \n", time);
		        nc = 0;
		        for (int l = 0; l < cvode_ncells ; ++l){
		            for (int sp=0;sp<Ncomp; sp++){
		                rhoY(indx_i[l],indx_j[l],indx_k[l],sp) = tmp_vect[l*(Ncomp+1) + sp];
		            }
		            rhoY(indx_i[l],indx_j[l],indx_k[l],Ncomp)  = tmp_vect[l*(Ncomp+1) + Ncomp];
		            rhoE(indx_i[l],indx_j[l],indx_k[l],0)      = tmp_vect_energy[l];
		        }
		    }
		}
	    }
	}
	if (nc != 0) {
		printf(" WARNING !! Not enough cells (%d) to fill %d \n", nc, cvode_ncells);
	} else {
		printf(" Integrated %d cells (4096)\n",num_cell_cvode_int);
	}
    }


    outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
    // Specs
    PlotFileFromMF(mf,outfile);
    
    extern_close();

    }

    amrex::Finalize();

    return 0;
}
